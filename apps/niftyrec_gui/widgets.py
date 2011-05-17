import pygtk
import gtk
from gtk import DrawingArea
import cairo

from StringIO import StringIO

def image2pixbuf(im):  
    file1 = StringIO()  
    im.save(file1, "ppm")  
    contents = file1.getvalue()  
    file1.close()  
    loader = gtk.gdk.PixbufLoader("pnm")  
    loader.write(contents, len(contents))  
    pixbuf = loader.get_pixbuf()  
    loader.close()  
    return pixbuf 


class ResizableImage(DrawingArea):
        def __init__(self, aspect_ratio=1.0, enlarge=False,
                interp=gtk.gdk.INTERP_HYPER, backcolor=None, max=(1600,1200), min=(20,15)):
            """Construct a ResizableImage control.
            Parameters:
            aspect_ratio -- Image sapect ratio (W/H) is multiplied by this factor. If 1.0 image keeps aspect ratio
            enlarge -- Allow image to be scaled up?
            interp -- Method of interpolation to be used.
            backcolor -- Tuple (R, G, B) with values ranging from 0 to 1,
                or None for transparent.
            max -- Max dimensions for internal image (width, height).
            """
            DrawingArea.__init__(self)
            self.pixbuf = None
            self.max = max
            self.min = min
            self.backcolor = backcolor
            self.interp = interp
            self.aspect_ratio = aspect_ratio
            self.enlarge = enlarge
            
            self.connect('expose_event', self.expose)
            self.add_events(gtk.gdk.KEY_PRESS_MASK | gtk.gdk.POINTER_MOTION_MASK | gtk.gdk.BUTTON_PRESS_MASK | gtk.gdk.SCROLL_MASK)

            self.x_pointer = 0
            self.y_pointer = 0
            self.color_pointer = (0.0,1.0,0.0)
            self._show_pointer = False
            self._width  = 0
            self._height = 0
            self._offset_w = 0
            self._offset_h = 0

        def set_pointer(self,x,y):
            self.x_pointer = (1.0*x-self._offset_w)/self._width
            self.y_pointer = (1.0*y-self._offset_h)/self._height
            self.invalidate()

        def set_pointer_01(self,x,y):
            self.x_pointer = x
            self.y_pointer = y
            self.invalidate()          

        def show_pointer(self,enable):
            self._show_pointer = enable
            self.invalidate()

        def set_color_pointer(self,rgb):
            self.color_pointer = rgb
            self.invalidate()

        def pointer_to_01(self,x,y):
            return ((1.0*x-self._offset_w)/self._width, (1.0*y-self._offset_h)/self._height)



        def set_aspect_ratio(self,aspect_ratio):
            self.aspect_ratio = aspect_ratio
            
        def expose(self, widget, event):
            # Load Cairo drawing context.
            self.context = self.window.cairo_create()
            # Set a clip region.
            self.context.rectangle(
                event.area.x, event.area.y,
                event.area.width, event.area.height)
            self.context.clip()
            # Render image.
            self.draw(self.context)
            return False


        def draw(self, context):
            # Get dimensions.
            rect = self.get_allocation()
            x, y = rect.x, rect.y
            # Remove parent offset, if any.
            parent = self.get_parent()
            if parent:
                offset = parent.get_allocation()
                x -= offset.x
                y -= offset.y
            # Fill background color.
            if self.backcolor:
                context.rectangle(x, y, rect.width, rect.height)
                context.set_source_rgb(*self.backcolor)
                context.fill_preserve()
            # Check if there is an image.
            if not self.pixbuf:
                return
            if rect.width < self.min[0] or  rect.height < self.min[1]:
                return
            width, height = resizeToFit(
                (self.pixbuf.get_width(), self.pixbuf.get_height()),
                (rect.width, rect.height),
                self.aspect_ratio,
                self.enlarge)
            self._width  = width
            self._height = height
            self._offset_w = (rect.width - width) / 2
            self._offset_h = (rect.height - height) / 2
            x = x + self._offset_w
            y = y + self._offset_h
            context.set_source_pixbuf(self.pixbuf.scale_simple(width, height, self.interp), x, y)
            context.paint()
            # Draw pointer
            if self._show_pointer:
                x_p = self.x_pointer*width
                y_p = self.y_pointer*height
                if x_p>width:
                    x_p=width
                if y_p>height:
                    y_p=height
                if x_p<0:
                    x_p=0
                if y_p<0:
                    y_p=0
                context.set_line_width(1.25)
                context.set_source_rgb(self.color_pointer[0], self.color_pointer[1], self.color_pointer[2])
                context.move_to(x_p+self._offset_w,self._offset_h)
                context.rel_line_to(0, height)
                context.move_to(self._offset_w, y_p+self._offset_h)
                context.rel_line_to(width, 0)
                context.stroke()


        # Handle the expose-event by drawing
        def do_expose_event(self, event):
            # Create the cairo context
            cr = self.window.cairo_create()
            # Restrict Cairo to the exposed area; avoid extra work
            cr.rectangle(event.area.x, event.area.y, event.area.width, event.area.height)
            cr.clip()
            self.draw(cr)

        def set_from_pixbuf(self, pixbuf):
            width, height = pixbuf.get_width(), pixbuf.get_height()
            # Limit size of internal pixbuf to increase speed.
            if (width > self.max[0] or height > self.max[1]):
                width, height = resizeToFit((width, height), self.max)
                self.pixbuf = pixbuf.scale_simple(width, height,gtk.gdk.INTERP_BILINEAR)                
            elif (width < self.min[0] or height < self.min[1]):
                width, height = resizeToFit((width, height), self.min)
                self.pixbuf = pixbuf.scale_simple(width, height,gtk.gdk.INTERP_BILINEAR)
            else:
                self.pixbuf = pixbuf
            self.invalidate()

        def set_from_file(self, filename):
            self.set_from_pixbuf(gtk.gdk.pixbuf_new_from_file(filename))

        def invalidate(self):
            self.queue_draw()


def resizeToFit(image, frame, aspect_ratio=1.0, enlarge=False):
        """Resizes a rectangle to fit within another.

        Parameters:
        image -- A tuple of the original dimensions (width, height).
        frame -- A tuple of the target dimensions (width, height).
        aspect_ratio -- Multiply image aspect ratio by this factor
        enlarge -- Allow image to be scaled up?

        """
        #if aspect:
        return scaleToFit(image, frame, aspect_ratio, enlarge)
        #else:
        #    return stretchToFit(image, frame, enlarge)

def scaleToFit(image, frame, aspect_ratio, enlarge=False):
        image_width, image_height = image
        frame_width, frame_height = frame
        image_aspect = float(image_width) / image_height
        frame_aspect = float(frame_width) / frame_height
        # Determine maximum width/height (prevent up-scaling).
        if not enlarge:
            max_width = min(frame_width, image_width)
            max_height = min(frame_height, image_height)
        else:
            max_width = frame_width
            max_height = frame_height
        # Frame is wider than image.
        if frame_aspect > (image_aspect*aspect_ratio):
            height = max_height
            width = int(height * (image_aspect*aspect_ratio))
        # Frame is taller than image.
        else:
            width = max_width
            height = int(width / (image_aspect*aspect_ratio))
        return (width, height)

#def stretchToFit(image, frame, enlarge=False):
#        image_width, image_height = image
#        frame_width, frame_height = frame
#        # Stop image from being blown up.
#        if not enlarge:
#            width = min(frame_width, image_width)
#            height = min(frame_height, image_height)
#        else:
#            width = frame_width
#            height = frame_height
#        return (width, height)


