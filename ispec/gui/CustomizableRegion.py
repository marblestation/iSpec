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

class CustomizableRegion:
    min_width = 0.002 # nm
    mark = None

    def __init__(self, frame, element_type, left_edge, right_edge, mark_position=None, note_text=None):
        self.frame = frame
        # Line specific:
        self.press = None
        self.note = None
        self.mark = None
        self.mark_position = mark_position

        self.element_type = element_type
        if self.element_type == "continuum":
            self.axvspan = self.frame.axes.axvspan(left_edge, right_edge, facecolor='green', alpha=0.30)
            self.axvspan.zorder = 3
        elif self.element_type == "lines":
            self.axvspan = self.frame.axes.axvspan(left_edge, right_edge, facecolor='yellow', alpha=0.30)
            self.axvspan.zorder = 3
            self.mark = self.frame.axes.axvline(x = self.mark_position, linewidth=1, color='orange')

            if note_text != "":
                self.note = self.frame.axes.annotate(note_text, xy=(self.mark_position, 1),  xycoords=("data", 'axes fraction'),
                    xytext=(-10, 20), textcoords='offset points',
                    size=8,
                    bbox=dict(boxstyle="round", fc="0.8"),
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                    edgecolor='black'),
                    horizontalalignment='right', verticalalignment='top',
                    annotation_clip=True,
                    )
            else:
                self.note = None

        elif self.element_type == "segments":
            self.axvspan = self.frame.axes.axvspan(left_edge, right_edge, facecolor='grey', alpha=0.30)
            # Segments always in the background but above spectra
            self.axvspan.zorder = 2
        self.original_facecolor = self.axvspan.get_facecolor()
        self.original_edgecolor = self.axvspan.get_edgecolor()

        # Fit line properties, dictionary of spectrum (to vinculate to different spectrum):
        self.line_plot_id = {}
        self.line_model = {}
        self.line_extra = {}

    def connect(self):
        # Connect to all the events
        self.cid_press = self.axvspan.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.axvspan.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.axvspan.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)


    def update_size(self, event):
        # Update the position of the edge left or right
        xy = self.axvspan.get_xy()

        # If left button clicked, modify left edge
        if event.button == 1:
            # Check condition to allow modification
            if self.mark is not None:
                # Mark should be inside the region
                compatible_with_mark = event.xdata < self.mark.get_xdata()[0]
            else:
                # There is no mark, so any position is good
                compatible_with_mark = True

            # Do not modify if the region will become slimmer than...
            big_enough = xy[2,0] - event.xdata > self.min_width

            if big_enough and compatible_with_mark:
                xy[0,0] = event.xdata
                xy[1,0] = event.xdata
                xy[4,0] = event.xdata
                self.frame.status_message("Moving left edge to %.4f" % xy[0,0])
                self.frame.canvas.draw()
        elif event.button == 3:
            # Check condition to allow modification
            if self.mark is not None:
                # Mark should be inside the region
                compatible_with_mark = event.xdata > self.mark.get_xdata()[0]
            else:
                # There is no mark, so any position is good
                compatible_with_mark = True

            # Do not modify if the region will become slimmer than...
            big_enough = event.xdata - xy[0,0] > self.min_width

            if big_enough and compatible_with_mark:
                xy[2,0] = event.xdata
                xy[3,0] = event.xdata
                self.frame.status_message("Moving right edge to %.4f" % xy[2,0])
                self.frame.canvas.draw()

    def update_mark(self, event):
        # Position of the edge left or right
        xy = self.axvspan.get_xy()

        # If left button clicked, modify mark
        if event.button == 1 and self.mark is not None:
            x = self.mark.get_xdata()

            inside_region = (event.xdata > xy[0,0]) and (event.xdata < xy[2,0])
            if inside_region:
                x[0] = event.xdata
                x[1] = event.xdata
                self.mark.set_xdata(x)
                self.mark_position = self.mark.get_xdata()[0]
                if self.note is not None:
                    self.note.xy = (event.xdata, 1)
                self.frame.status_message("Moving mark to %.4f" % x[0])
                self.frame.canvas.draw()

    def update_mark_note(self, event):
        if self.note is not None:
            note_text = self.note.get_text()
        else:
            note_text = ""

        note_text = self.frame.ask_value('Note for the new line region:', 'Note', note_text)
        self.update_mark_note_text(note_text)
        self.frame.canvas.draw()

    def update_mark_note_text(self, note_text):
        if note_text is not None and note_text != "":
            # Set
            if self.note is None:
                # New
                self.note = self.frame.axes.annotate(note_text, xy=(self.mark_position, 1),  xycoords=("data", 'axes fraction'),
                    xytext=(-10, 20), textcoords='offset points',
                    size=8,
                    bbox=dict(boxstyle="round", fc="0.8"),
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                    edgecolor='black'),
                    horizontalalignment='right', verticalalignment='top',
                    annotation_clip=True,
                    )
            else:
                # Update
                self.note.set_text(note_text)
        elif note_text == "" and self.note is not None:
            # Remove
            self.note.set_visible(False)
            self.note = None

    def show_mark(self):
        if self.mark is None:
            self.mark = self.frame.axes.axvline(x = self.mark_position, linewidth=1, color='orange')

    def hide_mark(self):
        if self.mark is not None:
            self.mark.remove()
            self.mark = None


    def on_press(self, event):
        # Validate that the click is on this axis/object
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes is not None and event.inaxes.get_navigate_mode() is not None: return
        contains, attrd = self.axvspan.contains(event)
        if not contains: return
        # This element is the kind of element that is selected to be modified?
        if self.frame.elements != self.element_type: return
        # If the action is "create", this should be managed by the frame and not individual elements
        if self.frame.action == "Create" and not (self.frame.elements == "lines" and self.frame.subelements == "marks"): return

        if self.frame.operation_in_progress:
            return

        # When regions overlap two or more can receive the click event, so
        # let's use a lock to allow the modification of one of them
        if self.frame.lock.acquire(False):
            if self.frame.action == "Remove":
                # The real removal is in the "on_release" event
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    self.press = event.button, event.x, event.xdata
                else:
                    self.press = event.button, event.x, event.xdata
                # Do not free the lock until the user releases the click
            elif self.frame.action == "Modify":
                ## Modify region
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata

                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    if event.button == 1:
                        # Left button, modify mark position
                        self.update_mark(event)
                    elif event.button == 2:
                        pass
                        # Modification of the note is managed on_release
                        # because we cannot show modal windows asking information
                        # while on_press or the mouse stops working
                else:
                    self.update_size(event)
                self.frame.regions_changed(self.element_type)
                # Do not free the lock until the user releases the click
            elif self.frame.action == "Create" and self.frame.elements == "lines" and self.frame.subelements == "marks":
                ## Create note but handel it on_release
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata
            elif self.frame.action == "Stats":
                ## Statistics about the region
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata


    def disconnect_and_remove(self):
        #self.hide()
        #self.frame.canvas.draw()
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_press)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_release)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_motion)
        self.axvspan.remove()
        if self.mark is not None:
            self.mark.remove()
        if self.note is not None:
            self.note.remove()

    def hide(self):
        self.axvspan.set_visible(False)
        if self.mark is not None:
            self.mark.set_visible(False)
        if self.note is not None:
            self.note.set_visible(False)
        if self.line_plot_id is not None:
            for plot_id in self.line_plot_id.values():
                if plot_id is not None:
                    self.frame.axes.lines.remove(plot_id)


    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes is not None and event.inaxes.get_navigate_mode() is not None: return

        # Validate it is not in PAN or ZOOM mode
        if event.inaxes is not None and event.inaxes.get_navigate_mode() is not None: return

        # This element is the kind of element that is selected to be modified?
        if self.frame.elements != self.element_type: return
        # If the action is "create", this should be managed by the frame and not individual elements
        if self.frame.action == "Create" and not (self.frame.elements == "lines" and self.frame.subelements == "marks"): return

        if self.press is not None:
            # In modification mode, if it is the current selected widget
            self.press = None
            if self.frame.action == "Remove":
                self.frame.lock.release() # Release now because the next commands will
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    # Remove note
                    if self.note is not None:
                        self.note.set_visible(False)
                        self.note.remove()
                        self.note = None
                        self.frame.canvas.draw()
                else:
                    # Remove from the figure now (except for line marks)
                    #  (if we do it on_press, this event would not be trigered and the lock not released)
                    self.frame.region_widgets[self.element_type].remove(self)
                    self.frame.regions_changed(self.element_type)
                    self.frame.flash_status_message("Removed region from " + "%.4f" % self.axvspan.get_xy()[0,0] + " to " + "%.4f" % self.axvspan.get_xy()[2,0])
                    self.disconnect_and_remove()
                    self.hide()
                    self.frame.canvas.draw()
            else:
                if self.frame.action == "Modify":
                    # If it was a right click when elements is lines and subelements marks
                    # => Change note
                    if event.button == 3 and self.frame.elements == "lines" and self.frame.subelements == "marks":
                        # Right button, modify note
                        if self.frame.action == "Modify" or self.frame.action == "Create":
                            self.update_mark_note(event)
                        else:
                            # Note removal is managed on_press
                            pass
                elif self.frame.action == "Create" and self.frame.elements == "lines" and self.frame.subelements == "marks":
                    # Create a note (left or right click)
                    if (event.button == 1 or event.button == 3) and self.note is None:
                        self.update_mark_note(event)
                elif self.frame.action == "Stats" and self.frame.lock.locked():
                    self.frame.update_stats(self)

                self.axvspan.set_facecolor(self.original_facecolor)
                self.axvspan.set_edgecolor(self.original_edgecolor)
                self.axvspan.set_alpha(0.30)
                self.frame.status_message("")
                self.frame.lock.release()
                self.frame.canvas.draw()

    def on_motion(self, event):
        # Validate that the object has been pressed and the click is on this axis
        if self.press is None: return
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes.get_navigate_mode() is not None: return
        if self.frame.action == "Stats": return

#        button, x, xpress = self.press
        if self.frame.action == "Modify":
            if self.frame.subelements == "marks":
                self.update_mark(event)
            else:
                self.update_size(event)


    def get_note_text(self):
        note_text = ""
        if self.element_type == "lines" and self.note is not None:
            note_text = self.note.get_text()
        return note_text

    def get_wave_peak(self):
        if self.element_type == "lines":
            return self.mark.get_xdata()[0]
        else:
            return None

    def get_wave_base(self):
        return self.axvspan.get_xy()[0,0]

    def get_wave_top(self):
        return self.axvspan.get_xy()[2,0]


