#!/usr/bin/env python

import pygtk
pygtk.require('2.0')
import gtk
import gtk.glade
import gobject

from ipython_view import *
import pango
import platform
if platform.system()=="Windows":
        FONT = "Lucida Console 9"
else:
        FONT = "Luxi Mono 10"

from widgets import ResizableImage, image2pixbuf

#if platform.system() == "Linux":
if(0):
    try:
        from volumeRender import VolumeRender
    except Exception,e:
        print e
        HAS_VOLUME_RENDER=False
    else:
        HAS_VOLUME_RENDER = True
else:
    try:
        from volumeRenderProxy import VolumeRender
    except Exception,e:
        print e
        HAS_VOLUME_RENDER=False
    else:
        HAS_VOLUME_RENDER = True

VOLUME_RENDER_DETECTOR_X = 512
VOLUME_RENDER_DETECTOR_Y = 512

from Reconstruction import Reconstructor

from Image import fromarray, new
from ImageDraw import Draw
from numpy import single, exp, uint16, double, ones

gtk.gdk.threads_init()
from time import sleep
from thread import start_new_thread
import os

NORM = 256

SCREENSHOT_FILENAME = './screenshot%d.png'


from colormap import colormap_store


class MainWindow:
    def __init__(self):
        #Set the Glade file
        self.gladefile = "niftyrec.glade"
        path = os.path.split(sys.argv[0])[0]
        self.wTree = gtk.glade.XML(path+"//"+self.gladefile) 
        #Get the Main Window, and connect the "destroy" event
        self.window = self.wTree.get_widget("MainWindow")
        if (self.window):
                self.window.connect("destroy", gtk.main_quit)
        self.window.connect("delete_event", self.delete_event)
        self.window.connect("destroy", self.destroy)
        connections = {"on_togglebutton_volume_render_toggled":self.on_togglebutton_volume_render_toggled,
                       "on_button_load_sinogram_clicked":self.on_click_load_sinogram,
                       "on_button_load_ct_clicked":self.on_click_load_ct,
                       "on_button_load_mri_clicked":self.on_click_load_mri,
                       "on_button_load_activity_clicked":self.on_click_load_activity, 
                       "on_button_clear_activity_clicked":self.on_click_clear_activity,
                       "on_button_view_activity_clicked":self.on_click_view_activity, 
                       "on_button_view_ct_clicked":self.on_click_view_ct,
                       "on_button_view_mri_clicked":self.on_click_view_mri,
                       "on_hscale_scale_sinogram_value_changed":self.on_changed_hscale_scale_sinogram,
                       "on_hscale_scale_viewport_value_changed":self.on_changed_hscale_scale_viewport,
                       "on_hscale_offset_sinogram_value_changed":self.on_changed_hscale_offset_sinogram,
                       "on_hscale_offset_viewport_value_changed":self.on_changed_hscale_offset_viewport,
                       "on_togglebutton_overlay_toggled":self.on_toggled_overlay,

                       "on_button_step_osem_clicked":self.on_click_step_osem, 
                       "on_button_step_tv_clicked":self.on_click_step_tv, 
                       "on_button_step_je_clicked":self.on_click_step_je, 
                       "on_button_step_ccje_clicked":self.on_click_step_ccje,
                       "on_button_reconstruct_osem_clicked":self.on_click_reconstruct_osem, 
                       "on_button_reconstruct_tv_clicked":self.on_click_reconstruct_tv, 
                       "on_button_reconstruct_je_clicked":self.on_click_reconstruct_je, 
                       "on_button_reconstruct_ccje_clicked":self.on_click_reconstruct_ccje,  
                       "on_button_stop_osem_clicked":self.on_click_stop_osem, 
                       "on_button_stop_tv_clicked":self.on_click_stop_tv, 
                       "on_button_stop_je_clicked":self.on_click_stop_je, 
                       "on_button_stop_ccje_clicked":self.on_click_stop_ccje,

                       "on_spinbutton_activity_size_x_value_changed":self.on_spinbutton_activity_size_x_value_changed,
                       "on_spinbutton_activity_size_y_value_changed":self.on_spinbutton_activity_size_y_value_changed,
                       "on_spinbutton_activity_size_z_value_changed":self.on_spinbutton_activity_size_z_value_changed,
                       "on_spinbutton_psf_fwhm_near_value_changed":self.on_spinbutton_psf_fwhm_near_value_changed,
                       "on_spinbutton_psf_fwhm_far_value_changed":self.on_spinbutton_psf_fwhm_far_value_changed,
                       "on_spinbutton_theta_first_value_changed":self.on_spinbutton_theta_first_value_changed,
                       "on_spinbutton_theta_last_value_changed":self.on_spinbutton_theta_last_value_changed,
                       "on_checkbutton_use_gpu_toggled":self.on_checkbutton_use_gpu_toggled,

                       "on_button_volume_render_preset_scene_1_clicked":self.on_button_volume_render_preset_scene_1_clicked, 
                       "on_button_volume_render_preset_scene_2_clicked":self.on_button_volume_render_preset_scene_2_clicked, 
                       "on_button_volume_render_preset_scene_3_clicked":self.on_button_volume_render_preset_scene_3_clicked, 
                       "on_button_volume_render_dump_screenshot_clicked":self.on_button_volume_render_dump_screenshot_clicked,
                      }
        self.wTree.signal_autoconnect(connections)
       
        self.progressbar = self.wTree.get_widget("progressbar")
    
        self.initializeScriptArea()
        self.initializeDrawingAreas()
        self.initializeFileList()
        self.initializeButtonsPanel()
        self.initializeVolumeRenderer()
        self.initializeGuiInputs()
        self.initializeReconstructor()
        self.initializeColormaps()
        self.exportNamespace()

        self.window.show()
        self.initialize_variables()

    def initialize_variables(self):
        self.displayed_image = None
        self.intensity_norm = 1
        self.intensity_scale = 1
        self.intensity_offset = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.intensity_norm_sinogram = 1
        self.intensity_scale_sinogram = 1
        self.intensity_offset_sinogram = 0
        self.n = 0
        self.show_overlay = False

    def delete_event(self, widget, event, data=None):
        # If you return FALSE in the "delete_event" signal handler,
        # GTK will emit the "destroy" signal. Returning TRUE means
        # you don't want the window to be destroyed.
        # This is useful for popping up 'are you sure you want to quit?'
        # type dialogs.
        print "delete event occurred"

        # Change FALSE to TRUE and the main window will not be destroyed
        # with a "delete_event".
        return False

    def destroy(self, widget, data=None):
        print "destroy signal occurred"
        gtk.main_quit()

    def main(self):
        gtk.gdk.threads_enter()
        gtk.main()
        gtk.gdk.threads_leave()

    def initializeScriptArea(self):
        self.script_window = gtk.ScrolledWindow()
        self.script_window.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
        self.python_view = IPythonView()
        self.python_view.modify_font(pango.FontDescription(FONT))
        self.python_view.set_wrap_mode(gtk.WRAP_CHAR)
        self.python_view.show()
        self.script_window.add(self.python_view)
        self.script_window.show()
        self.box_script = self.wTree.get_widget("box_script_")
        self.box_script.pack_start(self.script_window)
        self.box_script.show()

    def initializeDrawingAreas(self):
        self.box_drawing_area_1 = self.wTree.get_widget("box_drawing_area_1")
        self.draw_area_1 = ResizableImage(enlarge=True)
        self.box_drawing_area_1.add(self.draw_area_1)
        self.box_drawing_area_1.show()
        self.draw_area_1.show()
        self.draw_area_1.set_from_pixbuf(image2pixbuf(new('RGB',(10,10))))
#        self.draw_area_1.connect("button-press-event", wakeup)
        self.draw_area_1.connect("button-press-event", self.on_button_press_drawing_area_1)
        self.draw_area_1.connect("motion-notify-event", self.on_motion_drawing_area_1)
        self.draw_area_1.connect("scroll-event", self.on_scroll_drawing_area_1)
        self.draw_area_1.show_pointer(False)
        self.draw_area_1.set_pointer_01(0.5,0.5)
        self.draw_area_1.set_color_pointer((0.0,1.0,0.0))

        self.box_drawing_area_2 = self.wTree.get_widget("box_drawing_area_2")
        self.draw_area_2 = ResizableImage(enlarge=True)
        self.box_drawing_area_2.add(self.draw_area_2)
        self.box_drawing_area_2.show()
        self.draw_area_2.show()
        self.draw_area_2.set_from_pixbuf(image2pixbuf(new('RGB',(10,10))))
        self.draw_area_2.connect("button-press-event", self.on_button_press_drawing_area_2)
        self.draw_area_2.connect("motion-notify-event", self.on_motion_drawing_area_2)
        self.draw_area_2.connect("scroll-event", self.on_scroll_drawing_area_2)
        self.draw_area_2.show_pointer(False)
        self.draw_area_2.set_pointer_01(0.5,0.5)
        self.draw_area_2.set_color_pointer((0.0,1.0,0.0))

        self.box_drawing_area_3 = self.wTree.get_widget("box_drawing_area_3")
        self.draw_area_3 = ResizableImage(enlarge=True)
        self.box_drawing_area_3.add(self.draw_area_3)
        self.box_drawing_area_3.show()
        self.draw_area_3.show()
        self.draw_area_3.set_from_pixbuf(image2pixbuf(new('RGB',(10,10))))
        self.draw_area_3.connect("button-press-event", self.on_button_press_drawing_area_3)
        self.draw_area_3.connect("motion-notify-event", self.on_motion_drawing_area_3)
        self.draw_area_3.connect("scroll-event", self.on_scroll_drawing_area_3)
        self.draw_area_3.show_pointer(False)
        self.draw_area_3.set_pointer_01(0.5,0.5)
        self.draw_area_3.set_color_pointer((0.0,1.0,0.0))

        self.box_drawing_area_4 = self.wTree.get_widget("box_drawing_area_4")
        self.draw_area_4 = ResizableImage(enlarge=True)
        self.box_drawing_area_4.add(self.draw_area_4)
        self.box_drawing_area_4.show()
        self.draw_area_4.show()
        self.draw_area_4.set_from_pixbuf(image2pixbuf(new('RGB',(10,10))))
        self.draw_area_4.connect("button-press-event", self.on_button_press_drawing_area_4)
        self.draw_area_4.connect("motion-notify-event", self.on_motion_drawing_area_4)
        self.draw_area_4.connect("scroll-event", self.on_scroll_drawing_area_4)
        self.draw_area_4.show_pointer(False)
        self.draw_area_4.set_pointer_01(0.5,0.5)
        self.draw_area_4.set_color_pointer((1.0,0.0,0.0))

        # Timer to limit update rate:
        self._do_update_viewport = False
        gobject.timeout_add(30, self.on_viewport_timer)

        # Timer to initialize reconstructor:
        gobject.timeout_add(1000, self.on_initialize_reconstructor_timer)

    def on_viewport_timer(self):
        self._do_update_viewport = True
        return True

    def do_update_viewport(self):
        if self._do_update_viewport:
            self._do_update_viewport = False
            return True
        else:
            return False

    def on_initialize_reconstructor_timer(self):
        if self.size_has_changed:
            self.initializeReconstructor()
            print "Reconstructor reinitialized"
        return True



    def initializeColormaps(self):
        self.combobox_colormap1 = gtk.combo_box_new_text()
        self.box_colormap1 = self.wTree.get_widget("box_colormap1")
        self.box_colormap1.add(self.combobox_colormap1)
        self.combobox_colormap1.get_model().clear()
        self.combobox_colormap1.append_text("Colormap Functional")
        for colormap in colormap_store.keys():
            self.combobox_colormap1.append_text(colormap)
        self.combobox_colormap1.set_active(0)
        self.combobox_colormap1.connect("changed", self.on_colormap1_changed)
        self.combobox_colormap1.show()

        self.combobox_colormap2 = gtk.combo_box_new_text()
        self.box_colormap2 = self.wTree.get_widget("box_colormap2")
        self.box_colormap2.add(self.combobox_colormap2)
        self.combobox_colormap2.get_model().clear()
        self.combobox_colormap2.append_text("Colormap Anatomical")
        for colormap in colormap_store.keys():
            self.combobox_colormap2.append_text(colormap)
        self.combobox_colormap2.set_active(0)
        self.combobox_colormap2.connect("changed", self.on_colormap2_changed)
        self.combobox_colormap2.show()


    def initializeGuiInputs(self):
        self.spinbutton_activity_size_x = self.wTree.get_widget("spinbutton_activity_size_x")
        self.spinbutton_activity_size_y = self.wTree.get_widget("spinbutton_activity_size_y")
        self.spinbutton_activity_size_z = self.wTree.get_widget("spinbutton_activity_size_z")
        self.spinbutton_psf_fwhm_near = self.wTree.get_widget("spinbutton_psf_fwhm_near")
        self.spinbutton_psf_fwhm_far = self.wTree.get_widget("spinbutton_psf_fwhm_far")
        self.spinbutton_theta_first = self.wTree.get_widget("spinbutton_theta_first")
        self.spinbutton_theta_last = self.wTree.get_widget("spinbutton_theta_last")
        self.checkbutton_use_gpu = self.wTree.get_widget("checkbutton_use_gpu")
        self.spinbutton_subset_order_osem = self.wTree.get_widget("spinbutton_subset_order_osem")
        self.spinbutton_iterations_osem = self.wTree.get_widget("spinbutton_iterations_osem")
        self.spinbutton_subset_order_tv = self.wTree.get_widget("spinbutton_subset_order_tv")
        self.spinbutton_iterations_tv = self.wTree.get_widget("spinbutton_iterations_tv")
        self.spinbutton_beta_tv = self.wTree.get_widget("spinbutton_beta_tv")
        self.spinbutton_subset_order_je = self.wTree.get_widget("spinbutton_subset_order_je")
        self.spinbutton_iterations_je = self.wTree.get_widget("spinbutton_iterations_je")
        self.spinbutton_subset_order_ccje = self.wTree.get_widget("spinbutton_subset_order_ccje")
        self.spinbutton_iterations_ccje = self.wTree.get_widget("spinbutton_iterations_ccje")

    def initializeFileList(self):
        pass

    def initializeButtonsPanel(self):
        pass

    def initializeVolumeRenderer(self):
        self.volume_renderer = None
        self.activity_has_changed = False
        self.mri_has_changed = False

    def initializeReconstructor(self):
        if hasattr(self,'reconstructor'):
            if self.reconstructor.is_reconstructing():
                return 
        volume_size = (int(self.spinbutton_activity_size_x.get_value()),int(self.spinbutton_activity_size_y.get_value()),int(self.spinbutton_activity_size_z.get_value()))
        self.reconstructor=Reconstructor(volume_size)
        self.reconstructor.set_callback_status(self.set_reconstructor_status)
        self.reconstructor.set_callback_updateactivity(self.set_activity)
        self.progressbar.set_text(" ")
        self.progressbar.set_fraction(0.0)
        self.size_has_changed = False


    def exportNamespace(self):
        #self.python_view.updateNamespace({'NI_draw_area_1': self.draw_area_1})
        #self.python_view.updateNamespace({'NI_draw_area_2': self.draw_area_2})
        #self.python_view.updateNamespace({'NI_draw_area_3': self.draw_area_3})
        #self.python_view.updateNamespace({'NI_draw_area_4': self.draw_area_4})

        self.python_view.updateNamespace({'NI_set_view_1_from_image': self.set_view_1})
        self.python_view.updateNamespace({'NI_set_view_2_from_image': self.set_view_2})
        self.python_view.updateNamespace({'NI_set_view_3_from_image': self.set_view_3})
        self.python_view.updateNamespace({'NI_set_view_4_from_image': self.set_view_4})

        self.python_view.updateNamespace({'NI_set_coords_image_space': self.update_xyz_image})
        self.python_view.updateNamespace({'NI_set_coords_viewport_space': self.update_xyz})

        self.python_view.updateNamespace({'NI_volume_renderer': self.volume_renderer})
        self.python_view.updateNamespace({'NI_reconstructor': self.reconstructor})

    def on_togglebutton_volume_render_toggled(self, item, data=None):
        if not HAS_VOLUME_RENDER:
            return 
        if item.get_active():
            print "Volume Render"
            if self.volume_renderer != None:
                return
            else:
                volume_size = (int(self.spinbutton_activity_size_x.get_value()),int(self.spinbutton_activity_size_y.get_value()),int(self.spinbutton_activity_size_z.get_value()))
#                if platform.system()=='Linu':
                if(0):
                    self.volume_renderer = VolumeRender(volume_size,(VOLUME_RENDER_DETECTOR_X, VOLUME_RENDER_DETECTOR_Y))
                    self.volume_renderer.show()
                else:
                    print "Stating volumeRender process.." 
                    os.popen2('python volumeRender.py'+' '+str(volume_size[0])+' '+str(volume_size[1])+' '+str(volume_size[2])+' '+str(VOLUME_RENDER_DETECTOR_X)+' '+str(VOLUME_RENDER_DETECTOR_Y))
                    print "Creating volumeRender proxy.."
                    self.volume_renderer = VolumeRender()
                sleep(0.1)
                print "Updating volumeRender.."
 #               self.update_volume_renderer()
        else:
            if self.volume_renderer != None:
                self.volume_renderer.stop()
                self.volume_renderer = None
        self.python_view.updateNamespace({'NI_volume_renderer': self.volume_renderer})
    
    def dialog_get_image_filename(self,extensions,title):
        print "Browse DICOM integral"
        extensions = ["dcm","nii","nii.gz"]
        chooser = gtk.FileChooserDialog(title=title,action=gtk.FILE_CHOOSER_ACTION_OPEN,
                buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))
        filter = gtk.FileFilter();
        name = ''
        for extension in extensions:
            name = name + extension+" "
            filter.add_pattern("*."+extension)
        filter.set_name(name); 
        chooser.add_filter(filter)
        response = chooser.run()
        filename = ""
        if response == gtk.RESPONSE_OK:
            filename = str(chooser.get_filename())
            print filename, 'File Selected'                   
        elif response == gtk.RESPONSE_CANCEL:
            print 'No files selected' 
        chooser.destroy()
        return filename


    def on_click_load_sinogram(self,item,data=None):
        filename = self.dialog_get_image_filename(["dcm","nii","nii.gz"],"Load Sinogram")
        self.load_sinogram_from_file(filename)

    def on_click_load_ct(self,item,data=None):
        filename = self.dialog_get_image_filename(["dcm","nii","nii.gz"],"Load CT Image")
        self.load_ct_from_file(filename)

    def on_click_load_mri(self,item,data=None):
        filename = self.dialog_get_image_filename(["dcm","nii","nii.gz"],"Load MR Image")
        self.load_mri_from_file(filename)

    def on_click_load_activity(self,item,data=None):
        filename = self.dialog_get_image_filename(["dcm","nii","nii.gz"],"Load Activity Image")
        self.load_activity_from_file(filename)

    def on_click_clear_activity(self,item,data=None):
        self.reconstructor.clear_activity()
        self.set_activity(self.reconstructor.activity)

    def on_click_view_activity(self,item,data=None):
        self.set_display(self.reconstructor.activity,central_view=True)

    def on_click_view_ct(self,item,data=None):
        self.set_display(self.ct,central_view=True)

    def on_click_view_mri(self,item,data=None):
        self.set_display(self.mri,central_view=True)

    def on_changed_hscale_scale_sinogram(self,item,data=None):
        self.intensity_scale_sinogram = exp(item.get_value()/30)
        self.update_display_sinogram()

    def on_changed_hscale_scale_viewport(self,item,data=None):
        self.intensity_scale = exp(item.get_value()/30)
        self.update_display()

    def on_changed_hscale_offset_sinogram(self,item,data=None):
        self.intensity_offset_sinogram = NORM*(item.get_value()/100)
        self.update_display_sinogram()

    def on_changed_hscale_offset_viewport(self,item,data=None):
        self.intensity_offset = NORM*(item.get_value()/100)
        self.update_display()


    def on_button_press_drawing_area_1(self, item, data=None):
        (x,y) = self.draw_area_1.pointer_to_01(data.x,data.y)
        self.update_xyz(None,x,y)

    def on_motion_drawing_area_1(self,item,data=None):
        if (data.get_state() & gtk.gdk.BUTTON1_MASK):
            (x,y) = self.draw_area_1.pointer_to_01(data.x,data.y)
            self.update_xyz(None,x,y)

    def on_scroll_drawing_area_1(self,item,data=None):
        if data.direction == gtk.gdk.SCROLL_UP:
            x = self.x+1
            if x>=self.displayed_image.shape[0]-1:
                x = self.displayed_image.shape[0]
        else:
            x = self.x-1
            if x<0:
                x=0
        self.x=x
        self.update_display()

    def on_button_press_drawing_area_2(self, item, data=None):
        (x,y) = self.draw_area_2.pointer_to_01(data.x,data.y)
        self.update_xyz(x,None,y)

    def on_motion_drawing_area_2(self,item,data=None):
        if (data.get_state() & gtk.gdk.BUTTON1_MASK):
            (x,y) = self.draw_area_2.pointer_to_01(data.x,data.y)
            self.update_xyz(x,None,y)

    def on_scroll_drawing_area_2(self,item,data=None):
        if data.direction == gtk.gdk.SCROLL_UP:
            self.y+=1
            if self.y>=self.displayed_image.shape[1]-1:
                self.y = self.displayed_image.shape[1]
        else:
            self.y-=1
            if self.y<0:
                self.y=0
        self.update_display()

    def on_button_press_drawing_area_3(self, item, data=None):
        (x,y) = self.draw_area_3.pointer_to_01(data.x,data.y)
        self.update_xyz(y,x,None)

    def on_motion_drawing_area_3(self,item,data=None):
        if (data.get_state() & gtk.gdk.BUTTON1_MASK):
            (x,y) = self.draw_area_3.pointer_to_01(data.x,data.y)
            self.update_xyz(y,x,None)

    def on_scroll_drawing_area_3(self,item,data=None):
        if data.direction == gtk.gdk.SCROLL_UP:
            self.z+=1
            if self.z>=self.displayed_image.shape[2]-1:
                self.z = self.displayed_image.shape[2]
        else:
            self.z-=1
            if self.z<0:
                self.z=0
        self.update_display()

    def on_button_press_drawing_area_4(self, item, data=None):
        self.draw_area_4.set_pointer(data.x,data.y)

    def on_motion_drawing_area_4(self,item,data=None):
        if (data.get_state() & gtk.gdk.BUTTON1_MASK):
            self.draw_area_4.set_pointer(data.x,data.y)

    def on_scroll_drawing_area_4(self,item,data=None):
        if data.direction == gtk.gdk.SCROLL_UP:
            self.n+=1
            n_projections = self.sinogram.shape[0]
            if self.n>=n_projections-1:
                self.n=n_projections-1
        else:
            self.n-=1
            if self.n<0:
                self.n=0
        self.update_display_sinogram()

    def on_toggled_overlay(self, item ,data=None):
        self.enable_overlay(item.get_active())




    def on_click_step_osem(self, item, data=None):
        parameters = {'steps':1, 'subset_order':int(self.spinbutton_subset_order_osem.get_value())}
        self.reconstruct('osem',parameters)

    def on_click_reconstruct_osem(self, item, data=None):
        parameters = {'steps':int(self.spinbutton_iterations_osem.get_value()), 'subset_order':int(self.spinbutton_subset_order_osem.get_value())}
        self.reconstruct('osem',parameters)

    def on_click_stop_osem(self, item, data=None):
        self.reconstructor.stop()


    def on_click_step_tv(self, item, data=None):
        parameters = {'steps':1, 'beta':int(self.spinbutton_beta_tv.get_value()), 'subset_order':int(self.spinbutton_subset_order_tv.get_value())}
        self.reconstruct('tv',parameters)

    def on_click_reconstruct_tv(self, item, data=None):
        parameters = {'steps':int(self.spinbutton_iterations_tv.get_value()), 'beta':int(self.spinbutton_beta_tv.get_value()), 'subset_order':int(self.spinbutton_subset_order_tv.get_value())}
        self.reconstruct('tv',parameters)

    def on_click_stop_tv(self, item, data=None):
        self.reconstructor.stop()


    def on_click_step_je(self, item, data=None):
        parameters = {'steps':1, 'subset_order':int(self.spinbutton_subset_order_je.get_value())}
        self.reconstruct('je',parameters)

    def on_click_reconstruct_je(self, item, data=None):
        parameters = {'steps':int(self.spinbutton_iterations_je.get_value()), 'subset_order':int(self.spinbutton_subset_order_je.get_value())}
        self.reconstruct('je',parameters)

    def on_click_stop_je(self, item, data=None):
        self.reconstructor.stop()

    def on_click_step_ccje(self, item, data=None):
        parameters = {'steps':1, 'subset_order':int(self.spinbutton_subset_order_ccje.get_value())}
        self.reconstruct('ccje',parameters)

    def on_click_reconstruct_ccje(self, item, data=None):
        parameters = {'steps':int(self.spinbutton_iterations_ccje.get_value()), 'subset_order':int(self.spinbutton_subset_order_ccje.get_value())}
        self.reconstruct('ccje',parameters)

    def on_click_stop_ccje(self, item, data=None):
        self.reconstructor.stop()



    def on_spinbutton_activity_size_x_value_changed(self,item,data=None):
        x = item.get_value()
        self.spinbutton_activity_size_y.set_value(x)
        self.spinbutton_activity_size_z.set_value(x)
#        self.initializeReconstructor()
        self.size_has_changed = True

    def on_spinbutton_activity_size_y_value_changed(self,item,data=None):
        y = item.get_value()
        self.spinbutton_activity_size_x.set_value(y)
        self.spinbutton_activity_size_z.set_value(y)        
#        self.initializeReconstructor()
        self.size_has_changed = True

    def on_spinbutton_activity_size_z_value_changed(self,item,data=None):
        z = item.get_value()
        self.spinbutton_activity_size_x.set_value(z)
        self.spinbutton_activity_size_y.set_value(z)
#        self.initializeReconstructor()
        self.size_has_changed = True

    def on_spinbutton_psf_fwhm_near_value_changed(self,item,data=None):
        pass

    def on_spinbutton_psf_fwhm_far_value_changed(self,item,data=None):
        pass

    def on_spinbutton_theta_first_value_changed(self,item,data=None):
        pass

    def on_spinbutton_theta_last_value_changed(self,item,data=None):
        pass

    def on_checkbutton_use_gpu_toggled(self,item,data=None):
        pass


    def on_button_volume_render_preset_scene_1_clicked(self,item,data=None):
        pass

    def on_button_volume_render_preset_scene_2_clicked(self,item,data=None):
        pass

    def on_button_volume_render_preset_scene_3_clicked(self,item,data=None): 
        pass

    def on_button_volume_render_dump_screenshot_clicked(self,item,data=None):
        if not hasattr(self,"screenshot_number"): 
            self.screenshot_number = 0
        filename = SCREENSHOT_FILENAME%self.screenshot_number
        print "Saving image file: ",filename
        self.volume_renderer.dump_screenshot(filename)
        self.screenshot_number+=1


    def on_colormap1_changed(self,item,data=None):
        text = item.get_active_text()
        if text == "Colormap Functional": 
            return
        self.volume_renderer.set_colormap1(colormap_store[text])
 
    def on_colormap2_changed(self,item,data=None):
        text = item.get_active_text()
        if text == "Colormap Anatomical": 
            return
        self.volume_renderer.set_colormap2(colormap_store[text])



    def reconstruct(self,method,parameters):
        self.set_activity(self.reconstructor.activity,display_centered=True)
        N_cameras = self.sinogram.shape[0]
        self.reconstructor.set_cameras_equispaced(self.spinbutton_theta_first.get_value(),self.spinbutton_theta_last.get_value(),N_cameras,1)
        print "Use gpu: ", self.checkbutton_use_gpu.get_active()
        self.reconstructor.set_use_gpu(self.checkbutton_use_gpu.get_active())
        self.reconstructor.set_psf_two_points(self.spinbutton_psf_fwhm_near.get_value(), 0, self.spinbutton_psf_fwhm_far.get_value(),self.spinbutton_activity_size_z.get_value()-1)
        self.reconstructor.set_sinogram(self.sinogram)
#        self.reconstructor.set_sinogram(single(ones((128,128,120))))
#        self.reconstructor.set_attenuation(self.ct)
        self.reconstructor.reconstruct(method,parameters)

    def set_reconstructor_status(self,text,progress):
        self.progressbar.set_text(text)
        self.progressbar.set_fraction(progress)



    def set_display(self,image3D,central_view=False):
        self.displayed_image = image3D
        #set intensity scale
        self.intensity_norm = NORM/self.displayed_image.max()
        if central_view:
            s = self.displayed_image.shape
            self.x = s[0]/2
            self.y = s[0]/2
            self.z = s[0]/2
        self.update_display()

    def update_display(self):
        self.update_pointer()
        if not self.do_update_viewport():
            return False
        plane_x = fromarray(self.intensity_scale*(self.intensity_norm*self.displayed_image[self.x,:,:]-self.intensity_offset)).convert('RGB')
        plane_y = fromarray(self.intensity_scale*(self.intensity_norm*self.displayed_image[:,self.y,:]-self.intensity_offset)).convert('RGB')
        plane_z = fromarray(self.intensity_scale*(self.intensity_norm*self.displayed_image[:,:,self.z]-self.intensity_offset)).convert('RGB')
        #self.draw_area_1.set_from_pixbuf(gtk.gdk.pixbuf_new_from_array(plane_x))
        #self.draw_area_2.set_from_pixbuf(gtk.gdk.pixbuf_new_from_array(plane_y))
        #self.draw_area_3.set_from_pixbuf(gtk.gdk.pixbuf_new_from_array(plane_z))
        self.set_view_1(plane_x)
        self.set_view_2(plane_y)
        self.set_view_3(plane_z)

    def update_xyz(self,x_p,y_p,z_p):
        x=None; y=None; z=None
        if x_p!=None:
            x = int(x_p*self.displayed_image.shape[0])
        if y_p!=None:
            y = int(y_p*self.displayed_image.shape[1])
        if z_p!=None:
            z = int(z_p*self.displayed_image.shape[2])
        self.update_xyz_image(x,y,z)

    def update_xyz_image(self,x,y,z):
        if x!=None:
            self.x = int(x)
        if y!=None:
            self.y = int(y)
        if z!=None:
            self.z = int(z)
        self.update_display()
        
    def update_pointer(self):
        x_p = (1.0*self.x)/self.displayed_image.shape[0]
        y_p = (1.0*self.y)/self.displayed_image.shape[1]
        z_p = (1.0*self.z)/self.displayed_image.shape[2]
        self.draw_area_1.set_pointer_01(y_p,z_p)
        self.draw_area_2.set_pointer_01(x_p,z_p)
        self.draw_area_3.set_pointer_01(y_p,x_p)
 
    def update_display_sinogram(self):
        image = fromarray(self.intensity_scale_sinogram*(self.intensity_norm_sinogram*self.sinogram[self.n,:,:]-self.intensity_offset_sinogram)).convert('RGB')
        self.set_view_4(image)

    def update_volume_renderer(self):
        if not self.volume_renderer:
            return
        if self.activity_has_changed:
            self.volume_renderer.set_volume1(uint16(self.reconstructor.activity*0.9*2**16/self.reconstructor.activity.max()))
#        if self.mri_has_changed:
#            self.volume_renderer.set_volume2(self.mri)
        self.activity_has_changed=False
        self.mri_has_changed=False

    def load_image(self,filename):
        if filename.endswith('.dcm'):
            from dicom import read_file
            image = read_file(filename)
            image = single(image.pixel_array)
        elif (filename.endswith('.nii') or filename.endswith('.nii.gz')):
            from nifti import NiftiImage
            image = NiftiImage(filename)
            image = single(image.data)
        return image

    def load_sinogram_from_file(self,filename):
        image = self.load_image(filename)
        self.set_sinogram(image,display_centered=True)

    def load_activity_from_file(self,filename):
        image = self.load_image(filename)
        self.set_activity(image,display_centered=True)

    def load_ct_from_file(self,filename):
        image = self.load_image(filename)
        self.set_ct(image,display_centered=True)

    def load_mri_from_file(self,filename):
        image = self.load_image(filename)
        self.set_mri(image,display_centered=True)



    def set_sinogram(self,sinogram,display_centered=False):
        self.sinogram = sinogram
        self.n = self.sinogram.shape[0]/2
        self.intensity_norm_sinogram = NORM/self.sinogram.max()
        self.update_display_sinogram()
        self.python_view.updateNamespace({'NI_sinogram': self.sinogram})

    def set_activity(self,activity,display_centered=False):
        self.reconstructor.set_activity(activity)
        self.set_display(self.reconstructor.activity,display_centered)
        self.python_view.updateNamespace({'NI_activity': self.reconstructor.activity}) 
        self.python_view.updateNamespace({'NI_reconstructor': self.reconstructor}) 
        #update volume renderer
        self.activity_has_changed=True
        self.update_volume_renderer()

    def set_ct(self,ct,display_centered=False):
        self.ct = ct
        self.set_display(self.ct,display_centered)
        self.python_view.updateNamespace({'NI_ct': self.ct})

    def set_mri(self,mri,display_centered=False):
        self.mri_has_changed=True
        self.mri = mri
        self.set_display(self.mri,display_centered)
        self.python_view.updateNamespace({'NI_mri': self.mri})
        #update volume renderer
        self.mri_has_changed=True
        self.update_volume_renderer()



    def set_view_1(self,image):
        self.draw_area_1.set_from_pixbuf(image2pixbuf(image))

    def set_view_2(self,image):
        self.draw_area_2.set_from_pixbuf(image2pixbuf(image))

    def set_view_3(self,image):
        self.draw_area_3.set_from_pixbuf(image2pixbuf(image))

    def set_view_4(self,image):
        self.draw_area_4.set_from_pixbuf(image2pixbuf(image))


    def enable_overlay(self,enable):
        self.draw_area_1.show_pointer(enable)
        self.draw_area_2.show_pointer(enable)
        self.draw_area_3.show_pointer(enable)
        self.draw_area_4.show_pointer(enable)

# If the program is run directly or passed as an argument to the python
# interpreter then create a HelloWorld instance and show it
if __name__ == "__main__":
    hello = MainWindow()
    hello.main()

