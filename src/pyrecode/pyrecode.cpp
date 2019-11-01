/*#include <pyrecode_include.h>*/
#include <python.h>
#include "structmember.h"

#ifndef ENABLE_MULTIPLE_COMPRESSIONS
#define ENABLE_MULTIPLE_COMPRESSIONS 0
#endif


#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <omp.h>
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
#include "recode_utils.h"
#include "mrchandler.h"
#include "sequence_reader.h"
#include "commons.h"
#include "compressor.h"
#include "logger.h"
#include "fileutils.h"
#include "argparser.h"
#include "recode_header.h"

/*
#include "ReCoDeConfig.h"
*/

#include "L1.h"

typedef struct {
	PyObject_HEAD
	FILE *file;
	uint64_t current_frame_index;
} RecodeReader;

static PyMemberDef ReCoDe_members[] = {
	{ "file", T_OBJECT_EX, offsetof(RecodeReader, file), 0, "pointer to file" },
	{ "current_frame_index", T_INT, offsetof(RecodeReader, current_frame_index), 0, "index of next frame to be read" },
	{ NULL }  /* Sentinel */
};

static void 
ReCoDe_dealloc(RecodeReader *self)
{
	if (self->file != NULL) {
		fclose(self->file);
		Py_XDECREF(self->file);
	}
	Py_XDECREF(self->current_frame_index);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
ReCoDe_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	RecodeReader *self;
	self = (RecodeReader *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->file = NULL;
		self->current_frame_index = 0;
	}
	return (PyObject *)self;
}

PyObject * 
_open_file(RecodeReader* self, PyObject* args) {

	char * filename;
	if (!PyArg_ParseTuple(args, "s", &filename)) {
		return Py_BuildValue("i", 0);
	}
	self->file = fopen(filename, "r");
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

	// handle linux with fseeko

	uint64_t offset;
	int origin;
	if (!PyArg_ParseTuple(args, "kI", &offset, &origin)) {
		printf("0\n");
		return Py_BuildValue("i", 0);
	}
	printf("%d, %d\n", offset, origin);

	int state = _fseeki64 (self->file, offset, origin);
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

	/*To be moved to constructor*/
	uint32_t ny = 512;
	uint32_t nx = 4096;
	uint8_t bit_depth = 12;
	uint8_t byte_depth = ceil((bit_depth*1.0) / 8.0);

	uint32_t n_pixels_in_frame = nx * ny;														// number of pixels in the image
	uint32_t n_bytes_in_binary_image = ceil(n_pixels_in_frame / 8.0);							// number of bytes needed to pack binary image
	uint32_t n_bytes_in_image = ceil(n_pixels_in_frame * byte_depth);							// number of bytes needed to hold pixvals

	uint8_t *compressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));
	uint8_t *deCompressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));
	uint8_t *compressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	uint8_t *deCompressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));

	uint64_t pow2_lookup_table[64];
	uint8_t	 t;
	for (t = 0; t < 64; t++) {
		pow2_lookup_table[t] = pow(2, t);
	}
	/*To be moved to constructor*/

	uint32_t n_compressed_bytes_in_binary_image;
	uint32_t n_compressed_bytes_in_pixvals;
	uint32_t n_bytes_in_packed_pixvals;

	if (!PyArg_ParseTuple(args, "IIIy*", &n_compressed_bytes_in_binary_image, &n_compressed_bytes_in_pixvals, &n_bytes_in_packed_pixvals, &view_frameData)) {
		return Py_BuildValue("s", "Unable to parse argument frame_index");
	}
	printf("%d, %d, %d\n", n_compressed_bytes_in_binary_image, n_compressed_bytes_in_pixvals, n_bytes_in_packed_pixvals);

	decompressExpand_L1_Reduced_Compressed_Frame_Sparse (
		self->file, nx, ny, bit_depth, 
		n_compressed_bytes_in_binary_image, n_compressed_bytes_in_pixvals, n_bytes_in_packed_pixvals, n_bytes_in_binary_image, 
		compressedBinaryImage, deCompressedBinaryImage, compressedPixvals, deCompressedPixvals,
		pow2_lookup_table,
		(uint16_t *)(&view_frameData)->buf
	);
	return Py_BuildValue("i", 1);
}


PyMethodDef ReCoDeMethods[] = 
{
	{ "_compress_stream", (PyCFunction)_compress_stream, METH_VARARGS, 0 },
	{ "_decompress_stream_1", (PyCFunction)_decompress_stream_1, METH_VARARGS, 0 },
	{ "_open_file", (PyCFunction)_open_file, METH_VARARGS, 0 },
	{ "_close_file", (PyCFunction)_close_file, METH_VARARGS, 0 },
	{ "_fseek", (PyCFunction)_fseek, METH_VARARGS, 0 },
	{ "_get_frame_sparse_L1", (PyCFunction)_get_frame_sparse_L1, METH_VARARGS, 0 },
	{0,0,0,0}
};

static struct PyModuleDef PyReCoDe = {
	PyModuleDef_HEAD_INIT,
	"PyReCoDe",   /* name of module */
	"Interface to C ReCoDe Reader", /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
			  or -1 if the module keeps state in global variables. */
	ReCoDeMethods
};

static PyTypeObject ReCoDeReaderType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"PyReCoDe.Reader",                            /*tp_name*/
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
PyInit_PyReCoDe(void)
{
	//return PyModule_Create(&PyReCoDe);
	
	PyObject *m;
	if (PyType_Ready(&ReCoDeReaderType) < 0) {
		return NULL;
	}

	m = PyModule_Create(&PyReCoDe);
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