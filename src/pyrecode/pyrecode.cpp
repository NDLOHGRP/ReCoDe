#include <Python.h>
#include "structmember.h"

#ifndef ENABLE_MULTIPLE_COMPRESSIONS
#define ENABLE_MULTIPLE_COMPRESSIONS 0
#endif

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "zlib.h"

#if ENABLE_MULTIPLE_COMPRESSIONS
#include "bzlib.h"
#include "lzma.h"
#include "snappy-c.h"
#include "lz4.h"
#endif

#ifdef _WIN32
	#include <win32\dirent.h>
#else
	#include <dirent.h>
	#include <ctype.h>
#endif

#include "recodefs.h"
#include "compressor.h"

/*
#include "ReCoDeConfig.h"
*/

#include "reader.h"

typedef struct {
	PyObject_HEAD
	FILE *file;
	uint16_t ny;
	uint16_t nx;
	uint8_t bit_depth;
	uint8_t byte_depth;
	uint32_t n_pixels_in_frame;				// number of pixels in the image
	uint32_t n_bytes_in_binary_image;		// number of bytes needed to pack binary image
	uint32_t n_bytes_in_image;				// number of bytes needed to hold pixvals
	uint8_t *compressedBinaryImage;
	uint8_t *deCompressedBinaryImage;
	uint8_t *compressedPixvals;
	uint8_t *deCompressedPixvals;
	uint64_t pow2_lookup_table[64];
} RecodeReader;

static PyMemberDef ReCoDe_members[] = {
	{ "file", T_OBJECT_EX, offsetof(RecodeReader, file), 0, "pointer to file" },
	{ "ny", T_USHORT, offsetof(RecodeReader, ny), 0, "rows in frame" },
	{ "nx", T_USHORT, offsetof(RecodeReader, nx), 0, "cols in frame" },
	{ "bit_depth", T_UBYTE, offsetof(RecodeReader, bit_depth), 0, "bits per pixel" },
	{ "byte_depth", T_UBYTE, offsetof(RecodeReader, byte_depth), 0, "bytes per pixel" },
	{ "n_pixels_in_frame", T_UINT, offsetof(RecodeReader, n_pixels_in_frame), 0, "bytes per frame" },
	{ "n_bytes_in_binary_image", T_UINT, offsetof(RecodeReader, n_bytes_in_binary_image), 0, "bytes per frame" },
	{ "n_bytes_in_image", T_UINT, offsetof(RecodeReader, n_bytes_in_image), 0, "bytes per frame" },
	{ "compressedBinaryImage", T_OBJECT_EX, offsetof(RecodeReader, compressedBinaryImage), 0, "buffer for compressed binary image" },
	{ "deCompressedBinaryImage", T_OBJECT_EX, offsetof(RecodeReader, deCompressedBinaryImage), 0, "buffer for de-compressed binary image" },
	{ "compressedPixvals", T_OBJECT_EX, offsetof(RecodeReader, compressedPixvals), 0, "buffer for compressed pixel intensity values" },
	{ "deCompressedPixvals", T_OBJECT_EX, offsetof(RecodeReader, deCompressedPixvals), 0, "buffer for de-compressed pixel intensity values" },
	{ "pow2_lookup_table", T_ULONG, offsetof(RecodeReader, pow2_lookup_table), 0, "look-up table for bit unpacking" },
	{ NULL }  /* Sentinel */
};

static void 
ReCoDe_dealloc(RecodeReader *self)
{
	if (self->file != NULL) {
		fclose(self->file);
		Py_XDECREF(self->file);
	}
	/*
	Py_XDECREF(self->nx);
	Py_XDECREF(self->ny);
	Py_XDECREF(self->bit_depth);
	Py_XDECREF(self->byte_depth);
	Py_XDECREF(self->n_pixels_in_frame);
	Py_XDECREF(self->n_bytes_in_binary_image);
	Py_XDECREF(self->n_bytes_in_image);
	Py_XDECREF(self->compressedBinaryImage);
	Py_XDECREF(self->deCompressedBinaryImage);
	Py_XDECREF(self->compressedPixvals);
	Py_XDECREF(self->deCompressedPixvals);
	Py_XDECREF(self->pow2_lookup_table);
	*/
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
ReCoDe_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	RecodeReader *self;
	self = (RecodeReader *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->file = NULL;
		self->ny = 0;
		self->nx = 0;
		self->bit_depth = 0;
		self->byte_depth = 0;
		self->n_pixels_in_frame = 0;									// number of pixels in the image
		self->n_bytes_in_binary_image = 0;		// number of bytes needed to pack binary image
		self->n_bytes_in_image = 0;		// number of bytes needed to hold pixvals
		self->compressedBinaryImage = NULL;
		self->deCompressedBinaryImage = NULL;
		self->compressedPixvals = NULL;
		self->deCompressedPixvals = NULL;
		for (uint8_t t = 0; t < 64; t++) {
			self->pow2_lookup_table[t] = pow(2, t);
		}
	}
	return (PyObject *)self;
}

static PyObject *
_create_read_buffers (RecodeReader* self, PyObject* args) {

	uint16_t ny, nx;
	uint8_t bit_depth;
	if (!PyArg_ParseTuple(args, "HHb", &ny, &nx, &bit_depth)) {
		printf("Failed to parse arguments in RecodeReader constructor. Expected three: ny, nx, bit_depth\n");
		return Py_BuildValue("i", 0);
	}
	if (self != NULL) {
		self->ny = ny;
		self->nx = nx;
		self->bit_depth = bit_depth;
		self->byte_depth = ceil((self->bit_depth*1.0) / 8.0);
		self->n_pixels_in_frame = (uint32_t)self->nx * (uint32_t)self->ny;									// number of pixels in the image
		self->n_bytes_in_binary_image = ceil(self->n_pixels_in_frame / 8.0);		// number of bytes needed to pack binary image
		self->n_bytes_in_image = ceil(self->n_pixels_in_frame * self->byte_depth);		// number of bytes needed to hold pixvals
		self->compressedBinaryImage = (uint8_t*)calloc(self->n_bytes_in_binary_image, sizeof(uint8_t));
		self->deCompressedBinaryImage = (uint8_t*)calloc(self->n_bytes_in_binary_image, sizeof(uint8_t));
		self->compressedPixvals = (uint8_t*)calloc(self->n_bytes_in_image, sizeof(uint8_t));
		self->deCompressedPixvals = (uint8_t*)calloc(self->n_bytes_in_image, sizeof(uint8_t));
	}
	return Py_BuildValue("i", 1);
}

PyObject * 
_open_file(RecodeReader* self, PyObject* args) {

	char * filename;
	if (!PyArg_ParseTuple(args, "s", &filename)) {
		return Py_BuildValue("i", 0);
	}
	self->file = fopen(filename, "rb");
	if (self->file==NULL) {
		return Py_BuildValue("i", 0);
	}
	return Py_BuildValue("i", 1);
}

PyObject *
_close_file(RecodeReader* self) {
	int state = fclose(self->file);
	self->file = NULL;
	if (state==0) {
		return Py_BuildValue("i", 1);
	} else {
		return Py_BuildValue("i", 0);
	}
}

PyObject *
_fseek(RecodeReader* self, PyObject* args) {

	long offset;
	int origin;
	if (!PyArg_ParseTuple(args, "li", &offset, &origin)) {
		printf("0\n");
		return Py_BuildValue("i", 0);
	}
	printf("offset = %d, origin = %d\n", offset, origin);

	//Since we are seeking on rc file, fseek will do
	int state = fseek (self->file, offset, origin);
	//int state = _fseeki64 (self->file, 248073, origin);
	//int state = _fseeki64 (self->file, 248073, 0);
   	if (state)
    	return Py_BuildValue("i", 0);
	else
		return Py_BuildValue("i", 1);
}

/*
_decompress_stream (compression_scheme, compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen)
decompressed nCompressedLen bytes from buffer compressedData to buffer deCompressedData using compression_scheme and compression_level
the expected size of the decompressed data is nDataLen
*/
PyObject *
_compress_stream(RecodeReader *self, PyObject* args) {

	uint8_t  compression_scheme, compression_level;
	uint8_t  *data, *compressedData;
	uint32_t n_data_bytes, n_compressed_bytes;
	Py_buffer view_Data, view_CompressedData;

	if (!PyArg_ParseTuple(args, "BBy*Iy*", &compression_scheme, &compression_level, &view_Data, &n_data_bytes, &view_CompressedData)) {
		return Py_BuildValue("s", "Unable to parse argument frame_index");
	}
	compress_stream(compression_scheme, compression_level, (uint8_t *)(&view_Data)->buf, n_data_bytes, &n_compressed_bytes, (uint8_t *)(&view_CompressedData)->buf);
	printf("%d\n", n_compressed_bytes);

	return Py_BuildValue("i", n_compressed_bytes);
}

/*
_decompress_stream (compression_scheme, compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen)
decompressed nCompressedLen bytes from buffer compressedData to buffer deCompressedData using compression_scheme and compression_level
the expected size of the decompressed data is nDataLen
*/
PyObject *
_decompress_stream_1(RecodeReader *self, PyObject* args) {

	uint8_t  compression_scheme, compression_level;
	uint8_t  *compressedData, *deCompressedData;
	uint32_t nCompressedLen, nDataLen;
	Py_buffer view_compressedData, view_deCompressedData;

	if (!PyArg_ParseTuple(args, "BBy*y*II", &compression_scheme, &compression_level, &view_compressedData, &view_deCompressedData, &nCompressedLen, &nDataLen)) {
		return Py_BuildValue("s", "Unable to parse argument frame_index");
	}
	printf("%d, %d, %d, %d\n", compression_scheme, compression_level, nCompressedLen, nDataLen);
	decompress_stream(compression_scheme, compression_level, (uint8_t *)(&view_compressedData)->buf, (uint8_t *)(&view_deCompressedData)->buf, nCompressedLen, nDataLen);
	return Py_BuildValue("i", 1);
}

PyObject *
_get_frame_sparse_L1 (RecodeReader *self, PyObject* args) {

	Py_buffer view_frameData;
	uint32_t n_compressed_bytes_in_binary_image;
	uint32_t n_compressed_bytes_in_pixvals;
	uint32_t n_bytes_in_packed_pixvals;
	uint8_t is_intermediate;
	uint8_t get_frame_id;

	if (!PyArg_ParseTuple(args, "IIIbby*", &n_compressed_bytes_in_binary_image, &n_compressed_bytes_in_pixvals, &n_bytes_in_packed_pixvals, &is_intermediate, &get_frame_id, &view_frameData)) {
		return Py_BuildValue("s", "Unable to parse argument frame_index");
		return Py_BuildValue("k", 0);
	}
	/*
	printf("%d, %d, %d\n", n_compressed_bytes_in_binary_image, n_compressed_bytes_in_pixvals, n_bytes_in_packed_pixvals);
	printf("current position (pyrecode.cpp): %d\n", ftell(self->file));
	*/
	int64_t n = decompressExpand_L1_Reduced_Compressed_Frame_Sparse (
		self->file, self->nx, self->ny, self->bit_depth, 
		n_compressed_bytes_in_binary_image, n_compressed_bytes_in_pixvals, n_bytes_in_packed_pixvals, self->n_bytes_in_binary_image, 
		self->compressedBinaryImage, self->deCompressedBinaryImage, self->compressedPixvals, self->deCompressedPixvals,
		self->pow2_lookup_table,
		(uint16_t *)(&view_frameData)->buf,
		is_intermediate,
		get_frame_id
	);
	//printf("Decoded Frame with %d foreground pixels\n", n);
	return Py_BuildValue("L", n);
}


PyMethodDef ReCoDeMethods[] = 
{
	{ "_compress_stream", (PyCFunction)_compress_stream, METH_VARARGS, 0 },
	{ "_decompress_stream_1", (PyCFunction)_decompress_stream_1, METH_VARARGS, 0 },
	{ "_open_file", (PyCFunction)_open_file, METH_VARARGS, 0 },
	{ "_close_file", (PyCFunction)_close_file, METH_VARARGS, 0 },
	{ "_fseek", (PyCFunction)_fseek, METH_VARARGS, 0 },
	{ "_create_read_buffers", (PyCFunction)_create_read_buffers, METH_VARARGS, 0 },
	{ "_get_frame_sparse_L1", (PyCFunction)_get_frame_sparse_L1, METH_VARARGS, 0 },
	{0,0,0,0}
};

static struct PyModuleDef c_recode = {
	PyModuleDef_HEAD_INIT,
	"c_recode",   /* name of module */
	"Interface to C ReCoDe Reader", /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
			  or -1 if the module keeps state in global variables. */
	ReCoDeMethods
};

static PyTypeObject ReCoDeReaderType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"c_recode.Reader",                            /*tp_name*/
	sizeof(RecodeReader),                          		/*tp_basicsize*/
	0,											/*tp_itemsize*/
	(destructor)ReCoDe_dealloc,						/*tp_dealloc*/
	0,                                          /*tp_print*/
	0,                                          /*tp_getattr*/
	0,                                          /*tp_setattr*/
	0,                                          /*tp_compare*/
	0,                                          /*tp_repr*/
	0,                                          /*tp_as_number*/
	0,                                          /*tp_as_sequence*/
	0,                                          /*tp_as_mapping*/
	0,                                          /*tp_hash */
	0,                                          /*tp_call*/
	0,                                          /*tp_str*/
	0,                                          /*tp_getattro*/
	0,                                          /*tp_setattro*/
	0,                                          /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,								/*tp_flags*/
	"C Recode Reader object",                               /*tp_doc*/
	0,											/*tp_traverse*/
	0,											/*tp_clear*/
	0,                                          /*tp_richcompare*/
	0,                                          /*tp_weaklistoffset*/
	0,                                          /*tp_iter*/
	0,                                          /*tp_iternext*/
	ReCoDeMethods,                                 /*tp_methods*/
	ReCoDe_members,                                /*tp_members*/
	0,                                          /*tp_getsets*/
	0,                                          /*tp_base*/
	0,                                          /*tp_dict*/
	0,                                          /*tp_descr_get*/
	0,                                          /*tp_descr_set*/
	0,                                          /*tp_dictoffset*/
	0,											/*tp_init*/
	0,                                          /*tp_alloc*/
	ReCoDe_new,										/*tp_new*/
};


PyMODINIT_FUNC
PyInit_c_recode(void)
{
	//return PyModule_Create(&c_recode);
	
	PyObject *m;
	if (PyType_Ready(&ReCoDeReaderType) < 0) {
		return NULL;
	}

	m = PyModule_Create(&c_recode);
	if (m == NULL) {
		return NULL;
	}

	Py_INCREF(&ReCoDeReaderType);
	if (PyModule_AddObject(m, "Reader", (PyObject *)&ReCoDeReaderType) < 0) {
        Py_DECREF(&ReCoDeReaderType);
        Py_DECREF(m);
        return NULL;
    }

	return m;
	
}