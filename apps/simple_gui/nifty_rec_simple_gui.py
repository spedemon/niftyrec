#!/usr/bin/env python

glade_file = "nifty_rec_simple_gui.glade"
img_bg = "bg.png"

import sys
import os

try:
    import pygtk
    pygtk.require("2.0")
except:
    pass
try:
    import gtk
    import gtk.glade
except:
    sys.exit(1)
import gobject
import StringIO
from Image import open as open_image
from time import time, sleep
from NiftyRec import NiftyRec
gobject.threads_init()


glade_file = os.path.join(os.path.dirname(__file__), glade_file)
img_bg = os.path.join(os.path.dirname(__file__), img_bg)


def image2pixbuf(im,resize=False):
    if resize:
        im = im.resize(resize)
    file1 = StringIO.StringIO()
    im.save(file1,"ppm")
    contents = file1.getvalue()
    file1.close()
    loader = gtk.gdk.PixbufLoader("pnm")
    loader.write(contents,len(contents))
    pixbuf = loader.get_pixbuf()
    loader.close()
    return pixbuf



class MainWindow:
    """Main Window"""

    def __init__(self):
        #Set the Glade file
        self.gladefile = glade_file 
        self.wTree = gtk.glade.XML(self.gladefile)
	
        #Connect widgets
        dic = {"on_MainWindow_destroy" : gtk.main_quit,
               "on_spin_detector_changed": self._on_change_detector,
               "on_button_save_all_clicked": self._on_click_save_all,
               "on_button_compute_clicked":self._on_click_initialize_spect,
               "on_spin_n_detectors_changed": self._on_change_n_detectors,
               "on_button_backproject_clicked":self._on_backproject,
               "on_toggle_back_toggled":self._on_toggle_show_back_projection }
        self.wTree.signal_autoconnect(dic)

        self.main_window = self.wTree.get_widget("MainWindow")
        self.drawing_area = self.wTree.get_widget("image")
        self.spin_detector = self.wTree.get_widget("spin_detector")
        self.spin_n_detectors = self.wTree.get_widget("spin_n_detectors")
        self.progress_bar = self.wTree.get_widget("progress_bar")
        self.toggle_back = self.wTree.get_widget("toggle_back")
#        self.spin_detector.set_range(0,self.N_detectors-1+10)
#        self.spin_n_detectors.set_value(self.N_detectors)
        self.initializing = False
        self.initialized = False
        self.main_window.show()
        self.icon_draw()
        self.spect_imager = None
        self.inverse_previewer = None
        self.progress_bar.set_value(0)


    def icon_draw(self):
        image = open_image(img_bg)
        pixbuf = image2pixbuf(image,(3*image.size[0],3*image.size[1]))
        self.drawing_area.set_from_pixbuf(pixbuf) 
        pass

    def _on_click_save_all(self,widget,data=None):
        if self.initializing:
            return False
        if not self.initialized:
            return False
        #for i in range(self.N_detectors):
        #    image = self.spect_imager.get_projection(i,rotation_filter=BILINEAR)
        #    image = image_to_PIL(image).convert("RGB")   
        #    name = "./images/"+str(i)+".png"
        #    image.save(name)
        #    print name
        print "Done"   

    def _on_click_initialize_spect(self,widget,data=None):
        print "Initializing spect projector.."
        sleep(0.1)
        self.initializing = True
        #Store N detectors
#        self.N_detectors = int(self.spin_n_detectors.get_value())
#        self.spin_detector.set_range(0,self.N_detectors-1+10)
#        print "N detectors:",self.N_detectors
#        start_new_thread(self.initialize_spect,(self.N_detectors,))

    def _on_change_n_detectors(self,widget,data=None):
        print "N detectors changed to",int(widget.get_value())

    def _on_toggle_show_back_projection(self,widget,data=None):
        self._on_change_detector()

    def _on_change_detector(self,widget=None):
        print "Detector preview:",int(self.spin_detector.get_value())
        if self.initializing:
            return False
        if not self.initialized:
            return False
        self.change_detector(int(self.spin_detector.get_value()))

    def change_detector(self,index):
        if self.toggle_back.get_active():
        #if self.inverse_previewer:
            print "changing view:",
            return True
        #image = self.spect_imager.get_projection(index,rotation_filter=BILINEAR)
        #image = image_to_PIL(image).convert("RGB")
        #pixbuf = image2pixbuf(image,(3*image.size[0],3*image.size[1]))
        #self.drawing_area.set_from_pixbuf(pixbuf)            
        return True

    def _on_backproject(self,N):
        if not self.initialized:
            print "No projections"
            return False
        pass
        
    def backproject(self,N):
        self.progress_bar.set_value(1)
        t = time()
        self.progress_bar.set_value(0)
        return True

    def initialize_spect(self,N):
        self.progress_bar.set_value(3)
        pass



def main():
    mw = MainWindow()
    gtk.main()

if __name__ == "__main__":
    main()


