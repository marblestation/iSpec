#
#    This file is part of iSpec.
#    Copyright Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import os
import sampy
import tempfile
import signal
import sys
from astropy.io.votable import parse
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field, Group
import urllib.request, urllib.error, urllib.parse
import numpy as np
from astropy.io import fits as pyfits
import threading
import logging



class SAMPManager(object):
    def __init__(self, callback, check_connection_period=1):
        self.callback = callback
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.timer = None
        self.samp_hub = None
        self.samp_client = None
        self.metadata = {
                        'samp.name' : 'iSpec',
                        'samp.description.text' : 'iSpec',
                        'samp.icon.url' : 'file://'+self.dirname+'/images/iSpec.png',
                        }
        self.samp_client = sampy.SAMPIntegratedClient(metadata=self.metadata, addr='localhost')
        self.check_connection_period = check_connection_period # seconds
        self.__check_connection()
        signal.signal(signal.SIGINT, self.__signal_handler)

    def is_connected(self):
        working_connection = False
        if self.samp_client is not None and self.samp_client.isConnected():
            try:
                self.samp_client.ping()
                working_connection = True
            except:
                working_connection = False
        return working_connection

    def __check_connection(self):
        status = "Status: "
        if not self.is_connected():
            if self.samp_client is not None and self.samp_client.isConnected():
                # Destroy old client to completely reset the connection status
                del self.samp_client
                self.samp_client = sampy.SAMPIntegratedClient(metadata=self.metadata, addr='localhost')
                logging.info("SAMP Connection lost")
            status += "Disconnected"
            self.__connect()
        else:
            status += "Connected"
            pass
        logging.debug(status)

        self.timer = threading.Timer(self.check_connection_period, self.__check_connection)
        self.timer.start()

    def __connect(self):
        # Before we start, let's kill off any zombies
        try:
            self.shutdown()
        except:
            pass

        # sampy seems to fall over sometimes if 'localhost' isn't specified, even though it shouldn't
        #self.samp_hub = sampy.SAMPHubServer(addr='localhost')
        #self.samp_hub.start()

        try:
            self.samp_client.connect()
            logging.info("SAMP Connection established")
        except Exception:
            pass
        else:
            #samp_client.bindReceiveNotification("*", self.__samp_receive_notification)
            #samp_client.bindReceiveCall("*", self.__samp_receive_call)
            self.samp_client.bindReceiveNotification("samp.hub.event.register", self.__samp_receive_notification)
            self.samp_client.bindReceiveNotification("samp.hub.event.metadata", self.__samp_receive_notification)
            self.samp_client.bindReceiveNotification("samp.hub.event.subscriptions", self.__samp_receive_notification)
            self.samp_client.bindReceiveNotification("samp.hub.event.unregister", self.__samp_receive_notification)
            self.samp_client.bindReceiveNotification("table.load.votable", self.__samp_receive_notification)
            self.samp_client.bindReceiveNotification("spectrum.load.ssa-generic", self.__samp_receive_notification)
            self.samp_client.bindReceiveCall("table.load.votable", self.__samp_receive_call)
            self.samp_client.bindReceiveCall("spectrum.load.ssa-generic", self.__samp_receive_call)


    def __samp_receive_notification(self, private_key, sender_id, mtype, params, extra):
        #print "Notification:", private_key, sender_id, mtype, params, extra
        if mtype == 'samp.hub.event.register':
            #print "Registered:", params['id']
            pass
        elif mtype == 'samp.hub.event.metadata':
            #print "Metadata:", params['id'], "=", params['metadata']['samp.name']
            pass
        elif mtype == 'samp.hub.event.subscriptions':
            if 'spectrum.load.ssa-generic' in list(params['subscriptions'].keys()):
                #print params['id'], "supports spectrum"
                pass
            if 'table.load.votable' in list(params['subscriptions'].keys()):
                #print params['id'], "supports votable"
                pass
            #print params
        elif mtype == 'samp.hub.event.unregister':
            #print params['id'], "unregistered"
            pass
        else:
            # For instance, VOSpec sends load votable/spectrum in form of notification
            # so we should try to process them also
            msg_id = None
            spectrum = self.__samp_receive_and_transform_spectrum(mtype, params)
            logging.info("Spectrum received via SAMP")
            self.callback(spectrum, "Received_spectrum")

    # Function called when a call is received
    def __samp_receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
        #print "Call:", private_key, sender_id, msg_id, mtype, params, extra
        spectrum = self.__samp_receive_and_transform_spectrum(mtype, params)
        self.samp_client.ereply(msg_id, sampy.SAMP_STATUS_OK, result = {"txt": "printed"})
        logging.info("Spectrum received via SAMP")
        self.callback(spectrum, "Received_spectrum")

    def __samp_receive_and_transform_spectrum(self, mtype, params):
        if mtype == 'table.load.votable':
            #print "Received votable", params['url']
            votable = self.__read_votable(params['url'])
            spectrum =self. __votable_to_spectrum(votable)
            #print "Spectrum success"
        elif mtype == 'spectrum.load.ssa-generic':
            #print "Received spectrum", params['url'], "in format", params['meta']['Access.Format']
            if params['meta']['Access.Format'] == 'application/votable':
                #print "- Votable"
                votable = self.__read_votable(params['url'])
                spectrum = self.__votable_to_spectrum(votable)
                #print "Spectrum success"
            elif params['meta']['Access.Format'] == 'application/fits':
                #print "- FITS"
                spectrum = self.__read_vofits(params['url'])
                #print "Spectrum success"
            else:
                raise Exception("Unknown format")
        else:
            raise Exception("Unknown action")
        return spectrum

    def __signal_handler(self, signal, frame):
        print('SIGTERM received (ctrl+c)')
        self.shutdown()

    def shutdown(self):
        if self.timer is not None: self.timer.cancel()
        if self.samp_client.isConnected():
            try:
                self.samp_client.disconnect()
            except Exception:
                pass
        if self.samp_hub is not None and self.samp_hub._is_running: self.samp_hub.stop()
        sys.exit(0)

    def __get_target_client(self, target_client_name):
        neighbours = samp_client.getRegisteredClients()
        for neighbour in neighbours:
            metadata = self.samp_client.getMetadata(neighbour)
            try:
                if (metadata['samp.name'] == target_client_name):
                    return neighbour
            except KeyError:
                continue
        return None

    ########### [start] File operations
    def __votable_to_spectrum(self, votable):
        spectrum = None
        # Search for what it is required
        has_waveobs = False
        has_flux = False
        target_table = None
        waveobs_units = None
        for resource in votable.resources:
            for table in resource.tables:
                if len(table.array) > 0 and len(table.fields) >= 2:
                    for field in table.fields:
                        if field.name in ["wave", "waveobs"]:
                            has_waveobs = True
                            waveobs_units = field.unit
                        elif field.name == "flux":
                            has_flux = True
                        if has_waveobs and has_flux:
                            break
                if has_waveobs and has_flux:
                    target_table = table
                    break
            if has_waveobs and has_flux:
                break

        if target_table is not None:
            spectrum = np.recarray((len(target_table.array),), dtype=[('waveobs', float),('flux', float),('err', float)])
            if 'waveobs' in list(target_table.array.dtype.fields.keys()):
                spectrum['waveobs'] = target_table.array['waveobs']
            elif 'wave' in list(target_table.array.dtype.fields.keys()):
                spectrum['waveobs'] = target_table.array['wave']
            else:
                # the first column
                spectrum['waveobs'] = target_table.array[target_table.array.dtype.names[0]]
            if 'flux' in list(target_table.array.dtype.fields.keys()):
                spectrum['flux'] = target_table.array['flux']
            else:
                # the second column
                spectrum['flux'] = target_table.array[target_table.array.dtype.names[1]]
            if 'err' in list(target_table.array.dtype.fields.keys()):
                spectrum['err'] = target_table.array['err']
            elif 'error' in list(target_table.array.dtype.fields.keys()):
                spectrum['err'] = target_table.array['error']
            elif 'errors' in list(target_table.array.dtype.fields.keys()):
                spectrum['err'] = target_table.array['errors']
            if 'sigma' in list(target_table.array.dtype.fields.keys()):
                spectrum['err'] = target_table.array['sigma']
            elif len(table.fields) >= 3:
                # the third column if exists
                spectrum['err'] = target_table.array[target_table.array.dtype.names[2]]
            else:
                spectrum['err'] = 0.0
        if target_table is None or spectrum is None:
            raise Exception("Table not compatible")

        return spectrum

    def __read_votable(self, url):
        votable = None
        if url.startswith('http://'):
            u = urllib.request.urlopen(url)
            tmp = tempfile.NamedTemporaryFile(mode="wt", suffix=".xml", delete=False, encoding='utf-8')
            tmp.write(u.read())
            tmp.close()
            #print tmp.name
            votable = parse(tmp.name, pedantic=False)
            os.remove(tmp.name)
        elif url.startswith('file://localhost/'):
            filename = url[17:]
            if filename[0] != "/":
                filename = "/" + filename
            votable = parse(filename, pedantic=False)
        else:
            raise Exception("Unkown URL")
        return votable

    def __read_vofits(self, url):
        spectrum = None
        #print url
        if url.startswith('http://localhost/'):
            u = urllib.request.urlopen(url)
            tmp = tempfile.NamedTemporaryFile(mode="wt", suffix=".xml", delete=False, encoding='utf-8')
            tmp.write(u.read())
            tmp.close()
            #print tmp.name
            hdulist = pyfits.open(tmp.name)
            os.remove(tmp)
        elif url.startswith('file://localhost/'):
            filename = url[17:]
            if filename[0] != "/":
                filename = "/" + filename
            hdulist = pyfits.open(filename)
        else:
            raise Exception("Unkown URL")
        # Find hdu containing data
        target_hdu = None
        for hdu in hdulist:
            if hdu.data is not None and len(hdu.data) > 0 and len(list(hdu.data.dtype.fields.keys())) >= 2:
                target_hdu = hdu
        spectrum = np.recarray((len(target_hdu.data),), dtype=[('waveobs', float),('flux', float),('err', float)])
        # the first column
        spectrum['waveobs'] = target_hdu.data[target_hdu.data.dtype.names[0]]
        # the second column
        spectrum['flux'] = target_hdu.data[target_hdu.data.dtype.names[1]]
        if len(list(hdu.data.dtype.fields.keys())) >= 3:
            # the third column if exists
            spectrum['err'] = target_hdu.data[target_hdu.data.dtype.names[2]]
        else:
            spectrum['err'] = 0.0
        if target_hdu is None or spectrum is None:
            raise Exception("FITS not compatible")
        return spectrum

    def __spectrum_to_vofits(self, spectrum):
        t = pyfits.new_table(spectrum)
        t.header['TTYPE1'] = "WAVELENGTH"
        t.header['TTYPE2'] = "FLUX"
        t.header['TTYPE3'] = "SIGMA"
        fits = pyfits.HDUList(pyfits.PrimaryHDU())
        fits.append(t)
        return fits

    def __spectrum_to_votable(self, spectrum):
        # Create a new VOTable file...
        votable = VOTableFile()

        # ...with one resource...
        resource = Resource()
        votable.resources.append(resource)

        # ... with one table
        table = Table(votable)
        resource.tables.append(table)

        # Define some fields
        waveobs = Field(votable, name="waveobs", datatype="double", unit="nm", ucd="em.wl")
        flux = Field(votable, name="flux", datatype="double", unit="Jy", ucd="phot.flux")
        err = Field(votable, name="err", datatype="double", ucd="stat.error;phot.flux")
        table.fields.extend([waveobs, flux, err])
        table.groups.extend([Group([flux, err])])
        #import ipdb
        #ipdb.set_trace()
        # Now, use those field definitions to create the numpy record arrays, with
        # the given number of rows
        table.create_arrays(len(spectrum))

        # Now table.array can be filled with data
        table.array['waveobs'] = spectrum['waveobs']
        table.array['flux'] = spectrum['flux']
        table.array['err'] = spectrum['err']

        #votable.set_all_tables_format('binary') # VOSpec does not understand binary format
        return votable
    ########### [end]   File operations

    def get_subscribers(self):
        ids = []
        names = []
        as_tables = []

        try:
            if self.samp_client.isConnected():
                subscribers1 = self.samp_client.getSubscribedClients("table.load.votable")
                subscribers2 = self.samp_client.getSubscribedClients("spectrum.load.ssa-generic")
                subscribers = dict(list(subscribers1.items()) + list(subscribers2.items()))
                for subscriber in list(subscribers.keys()):
                    metadata = self.samp_client.getMetadata(subscriber)
                    name = ""
                    if 'samp.name' in metadata:
                        name = metadata['samp.name']
                    ids.append(subscriber)
                    names.append(name)
                    # Prefered form: "spectrum.load.ssa-generic"
                    if subscriber in list(subscribers2.keys()):
                        as_tables.append(False)
                    else:
                        as_tables.append(True)

                # If there are repeated names, add the identifier as a prefix
                uniq_names, uniq_names_index = np.unique(names, return_index=True)
                if len(names) != len(uniq_names):
                    for i in np.arange(len(names)):
                        # Only add a preffix to the repeated application names
                        if i not in uniq_names_index:
                            names[i] += " [" + ids[i] + "]"
        except Exception:
            pass

        return (ids, names, as_tables)

    # Broadcast a table file to TOPCAT
    def broadcast_spectrum(self, spectrum, spectrum_name, target_client, target_client_is_name = False, as_table=True, as_fits=False):
        if as_fits:
            suffix = ".fits"
        else:
            suffix = ".xml"
        tmp = tempfile.NamedTemporaryFile(mode="wt", suffix=suffix, delete=False, encoding='utf-8')
        tmp.close()

        if not as_fits:
            votable = self.__spectrum_to_votable(spectrum)
            votable.to_xml(tmp.name)
        else:
            fits = self.__spectrum_to_vofits(spectrum)
            fits.writeto(tmp.name, clobber=True)


        if as_table:
            mtype = 'table.load.votable'
        else:
            mtype = 'spectrum.load.ssa-generic'

        metadata = {
                    'samp.mtype' : mtype,
                    'samp.params' : {
                                     'name' : spectrum_name,
                                     'table-id' : spectrum_name,
                                     'url' : 'file://localhost/' + tmp.name
                                    }
                   }

        if target_client_is_name:
            target_client_name = target_client
            target_client = self.__get_target_client(target_client_name)
        response = False
        if target_client is not None:
            try:
                response = self.samp_client.notify(target_client, metadata)
                logging.info("Spectrum sent via SAMP")
            except Exception:
                pass
        # Do not delete or the target client will not have time to read it
        #os.remove(tmp.name)
        return response

#s = SAMPManager()
#ids, names, as_tables = s.get_subscribers()
#print names

#if len(names) > 0:
    #spectrum = ispec.read_spectrum("input/spectra/examples/narval_sun.s.gz")
    #print "*"
    #i = 0
    #s.broadcast_spectrum(spectrum, "test1", ids[i], as_table=as_tables[i])
    #print "1"
    ##j = 3
    ##s.broadcast_spectrum(spectrum, "test2", ids[j], as_table=as_tables[j])
    #print "end"

#signal.pause()
