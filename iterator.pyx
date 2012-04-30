cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void* malloc(size_t size)
    void free(void* mem)

from multiprocessing import Process, Queue, cpu_count

cdef class Iterator:
    cdef int W, H, max, strsize, cpus
    cdef double* real_axis
    cdef double* imag_axis
    cdef unsigned char* image_string
    cdef BL, TR, param, queue
    
    def __cinit__(self, int W=500, int H=500, int max=255,
                  BL=-1-1j, TR=1+1j, param=0+0j):
        if W <= 0 or H <= 0:
            self.image_string = self.real_axis = self.imag_axis = NULL
            return
        self.strsize = strsize = W*H
        self.image_string = <unsigned char *>malloc(strsize*sizeof(char))
        self.real_axis = <double *>malloc(W*sizeof(double));
        self.imag_axis = <double *>malloc(H*sizeof(double));

    def __dealloc__(self):
        free(self.image_string)
        free(self.real_axis)
        free(self.imag_axis)

    def __init__(self, int W=500, int H=500, int max=255,
                 BL=-1-1j, TR=1+1j, param=0+0j):
        if W <= 0 or H <= 0:
            raise ValueError("Image dimensions must be positive.")
        if max < 0 or max > 255:
            raise ValueError('The value of max must lie between 0 and 255.')
        if (self.image_string == NULL or
            self.real_axis == NULL or
            self.imag_axis == NULL):
            raise RuntimeError('Out of Memory')
        self.W, self.H, self.max = W, H, max
        self.cpus = cpu_count()
        self.queue = Queue()
        self.set(BL, TR, param)
        
    def __repr__(self):
        return 'Iterator: W=%s, max=%s, max=%s, BL=%s, TR=%s, param=%s'%(
            self.W, self.H, self.max, self.BL, self.TR, self.param)
    
    def build_axes(self):
        cdef double R, I, delta
        cdef int i
        R = self.BL.real
        delta = (self.TR.real - R)/self.W
        for i in range(self.W):
            self.real_axis[i] = R
            R += delta
        I = self.TR.imag # For graphics, the y axis points down
        delta = (self.BL.imag - I)/self.H
        for i in range(self.H):
            self.imag_axis[i] = I
            I += delta
        pass

    def set(self, BL, TR, param):
        self.BL, self.TR, self.param = BL, TR, param
        self.build_axes()

    def get_image(self):
        processes = []
        bands = list(range(self.cpus))
        result = bytes('')
        # Start a worker process for each cpu.
        for n in range(self.cpus):
            P = Process(target=self.iterate, args=[n, self.cpus])
            processes.append(P)
            P.start()
        # Collect the results
        for n in range(self.cpus):
            k, image = self.queue.get()
            bands[k] = image
        # Join all the subprocesses - "probably good practice"
        for P in processes:
            P.join()
        # Merge the bands 
        for band in bands:
            result += band
        return result

    def Xget_image(self):
        self.iterate(self.image_string, 0, 2)
        self.iterate(self.image_string, 1, 2)
        return bytes(self.image_string[:self.strsize])

    def iterate(self):
        # Subclasses override this method
        pass
    
cdef class Z_Iterator(Iterator):
    """
    Loops through all values in the specified band of the Z-rectangle,
    computing the number of steps before escape, using Z as the start.
    (Computes the Julia set for C = param.)
    """
    def iterate(self, int band=0, int num_bands=1):
        cdef int N, H0, H1, d, i, j, iterations, maxit=self.max+1
        cdef double R, I, RR, II
        cdef double Cr=self.param.real, Ci=self.param.imag
        cdef double *Zr = self.real_axis, *Zi=self.imag_axis
        cdef unsigned char* imgstr = self.image_string
        d = self.H/num_bands
        H0 = band*d
        H1 = self.H if band == num_bands-1 else H0 + d
        N = H0*self.W
        for j in range(H0, H1):
            for i in range(self.W):
                R, I = Zr[i], Zi[j]
                for iterations in range(maxit):
                    RR, II = R*R, I*I
                    if RR + II > 4.0:
                        break
                    I = 2*I*R + Ci
                    R = RR - II + Cr
                imgstr[N] = iterations
                N += 1
        self.queue.put( (band, bytes(imgstr[H0*self.W: N]) ) )
        
cdef class C_Iterator(Iterator):
    """
    Loops through all values in the specified band of the C-rectangle,
    computing the number of steps before the point specified as param
    escapes.  (Computes the Mandelbrot set, with param=0+0j.)
    """
    def iterate(self, int band=0, int num_bands=1):
        cdef int N, H0, H1, d, i, j, iterations, maxit=self.max+1
        cdef double R, I, RR, II
        cdef double R0=self.param.real, I0=self.param.imag
        cdef double *Cr=self.real_axis, *Ci=self.imag_axis
        cdef unsigned char* imgstr = self.image_string
        # Actually, each process gets its own copy of self.image_string
        # but I am having each one write to its own part, in case I
        # decide it is worthwhile trying to share the string instead
        # of using a queue.
        d = self.H/num_bands
        H0 = band*d
        H1 = self.H if band == num_bands-1 else H0 + d
        N = H0*self.W # may as well start here
        for j in range(H0, H1):
            for i in range(self.W):
                R, I = R0, I0
                for iterations in range(maxit):
                    RR, II = R*R, I*I
                    if RR + II > 4.0:
                        break
                    I = 2*I*R + Ci[j]
                    R = RR - II + Cr[i]
                imgstr[N] = iterations
                N += 1
        self.queue.put( (band, bytes(imgstr[H0*self.W: N]) ) )

def iterate(int W, int H, int max, z0, z1, c0, c1):
    """
    Generates an image (as a string) by iterating z -> z^2 + c.
    Arguments:
        W       image width
        H       image height
        max     maximum number of iterates (0 <= max <= 255)
        z0, z1  corners of a box in the z-domain
        c0, c1  corners of a box in the c-domain

    The boxes in the z- and c- domains are sampled on a grid of size
    WxH.  For each pair (i,j), the (i,j) grid point in the z-domain is
    used as an initial z-value and the (i,j) grid point in the
    c-domain is used as the value of c.  The function z -> z^2 + c is
    iterated up to max times, stopping if |z| > 2.  The (i,j) pixel in
    the image is set to the number of iterations performed.
    """
    cdef unsigned char *image_string
    cdef double *z_real, *z_imag, *c_real, *c_imag
    cdef double delta_z, delta_c, z, c
    cdef double R, I, RR, II, mod
    cdef int i, j, N, iterations, strsize=W*H

    # Validate arguments
    if W <= 0 or H <= 0:
        raise ValueError, "Image dimensions must be positive."
    if max < 0 or max > 255:
        raise ValueError('The value of max must lie between 0 and 255.')

    # Allocate memory
    image_string = <unsigned char *>malloc(strsize*sizeof(char))
    z_real = <double *>malloc(W*sizeof(double));
    z_imag = <double *>malloc(H*sizeof(double));
    c_real = <double *>malloc(W*sizeof(double));
    c_imag = <double *>malloc(H*sizeof(double));
    if (image_string == NULL or
        z_real == NULL or
        z_imag == NULL or
        c_real == NULL or
        c_imag == NULL):
        raise RuntimeError('Out of memory')

    # Initialize domains
    z, c = z0.real, c0.real
    delta_z, delta_c = (z1.real - z)/W, (c1.real - c)/W
    for i in range(W):
      z_real[i] = z
      z += delta_z
      c_real[i] = c
      c += delta_c
    z, c = z1.imag, c1.imag
    delta_z, delta_c = (z0.imag - z)/H, (c0.imag - c)/H
    for j in range(H):
      z_imag[j] = z
      z += delta_z
      c_imag[j] = c
      c += delta_c

    # Iterate
    N = 0
    for j in range(H):
        for i in range(W):
            I, R = z_imag[j], z_real[i]
            for iterations in range(max+1):
                RR, II = R*R, I*I
                if RR + II > 4.0:
                    break
                I = 2*I*R + c_imag[j]
                R = RR - II + c_real[i]
            image_string[N] = iterations
            N += 1

    result = bytes(image_string[:strsize])
    # Free memory
    free(z_real)
    free(z_imag)
    free(c_real)
    free(c_imag)
    free(image_string)
    # Return a string (as bytes)
    return result

def boundary(image, int W, int H, int max):
    """
    Generates a black and white image of the boundary of a colored region.
    Arguments:
        image    the image as a string
        W        the width of the image
        H        the height of the image
        max      threshold
    Returns an image which is the boundary of the region in which pixel
    values are less than max.
    """
    cdef unsigned char *in_string = image 
    cdef unsigned char *out_string
    cdef int length, i, j, p
    length = W*H
    # Validate arguments
    if length != len(image):
        raise ValueError('Invalid image dimensions.')
    if max < 0 or max > 255:
        raise ValueError('The value of max must lie between 0 and 255.')
    # Allocate memory
    out_string = <unsigned char *>malloc(length*sizeof(char))
    # Find the boundary
    for i in range(W):
        for j in range(H):
            p = i+j*W;
            out_string[p] = 255;
            if (in_string[p] < max):
                if ( (i > 0 and in_string[p-1] >= max) or
                     (i < W and in_string[p+1] >= max) or
                     (j > 0 and in_string[p-W] >= max) or
                     (j < H and in_string[p+W] >= max) ):
                    out_string[p] = 0;
    result = bytes(out_string[:length])
    # Free memory
    free(out_string)
    # Return a string (as bytes)
    return result

def boxcount(image, int W, int H, int max):
    """
    Inputs an image as a string. Counts how many boxes of size 2^i
    meet the image.  Returns a list of integers, where the ith integer
    is the number of boxes of size 2^i that meet the image.  The
    maximum box size is the largest power of 2 which is smaller than
    min(W,H)/4.  A box is deemed to meet the image if it contains a
    pixel of value >= max.

    Arguments:
        image    the image as a string
        W        the width of the image
        H        the height of the image
        max      threshold
    """
    cdef unsigned char *in_string = <unsigned char *>image
    cdef char *buf
    cdef int i, j, m, n, length, box_size
    cdef long count
    length = W*H
    # Validate arguments
    if length != len(image):
        raise ValueError('Invalid image dimensions.')
    if max < 0 or max > 255:
        raise ValueError('The value of max must lie between 0 and 255.')
    # Allocate memory
    buf = <char *>malloc(length*sizeof(char));
    # Compute dimension
    result = []
    for i in range(length):
        buf[i] = 1 if in_string[i] < max else 0
    box_size = (W if W < H else H) >> 2
    while box_size > 1:
        i, j = 0, 0
        for m in range(H, 0, -2):
            for n in range(W, 0, -2):
                if n > 1:
                    buf[i] = buf[j] | buf[j+1]
                    if m > 1:
                        buf[i] |= (buf[j+W] | buf[j+W+1])
                    j += 2
                else:
                    buf[i] = buf[j]
                    if m > 1:
                        buf[i] |= buf[j+W]
                    j += 1
                i += 1
            j += W;
        W = (W+1)/2
        H = (H+1)/2
        count = 0
        while i > 0:
            i -= 1
            count += buf[i]
        result.append(count)
        box_size = box_size >> 1;
    free(buf);
    return result;
