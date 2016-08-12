import vdif_assembler_cython
from multiprocessing import Process

class assembler:
	def __init__(self, arg1, arg2, flag, n):
		self._assembler = vdif_assembler_cython.assembler(arg1, arg2, flag, n)
		self.active = True

	def register_processor(self, p):
		self.python_processor = p

	def run(self):
		assembler_thread = Process(
	            target=self._assembler.run,
	          	)
		print "starting thread!"
		assembler_thread.start()
		while(True):
			print "getting chunk from assembler"
			t0, chunk = self._assembler.get_intensity_chunk()
			print "got chunk with t0: ", t0
			p.process_chunk(t0,1024, chunk, None)
		

class processor(object):
	"""
	To define a python processor, you subclass this base class.
	When the assembler runs, it will call process_chunk() with a sequence of chunks, represented
	by a (t0,nt,efield,mask) quadruple.
	Each chunk corresponds to range of timestamps [t0,t0+nt), where t0 is a 64-bit wraparound-free
	timestamp.
	WARNING 1: Usually these ranges will be contiguous between calls, e.g.
		[t0,t0+nt)   [t0+nt,t0+2*nt)   [t0+2*nt,t0+3*nt)   ...
	but the vdif_processor should not assume that this!  If there is a temporary 
	interruption in data stream, then a timestamp gap will appear.
	The 'efield' arg is a shape (nfreq,2,nt) complex array with electric field values, where
	the middle index is polarziation.  Missing data is represented by (0+0j).  The 'mask' arg
	is a shape (nfreq,2,nt) integer array which is 0 for missing data, and 1 for non-missing.
	 WARNING 2: Handling missing data is an important aspect of the vdif_processor since it 
	 happens all the time.  If a GPU correlator node is down, which is a frequent occurrence, 
	 then some frequencies will be "all missing".  There are also routine packet loss events 
	 on second-timescales which result in some high-speed samples being flagged as missing data.
	 """	
	def process_chunk(self, t0, nt, efield, mask):
		print 'process_chunk called! t0=%s nt=%s efield (%s,%s) mask (%s,%s)' % (t0, nt, efield.dtype, efield.shape, mask.dtype, mask.shape)



	def finalize(self):
		pass