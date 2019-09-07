/*#include <pyrecode_include.h>*/
#include <python.h>
#include "structmember.h"
#include <inttypes.h>
#include "zlib.h"
#include "compressor.h"

/*
_decompress_stream (compression_scheme, compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen)
decompressed nCompressedLen bytes from buffer compressedData to buffer deCompressedData using compression_scheme and compression_level
the expected size of the decompressed data is nDataLen
*/
PyObject *
_compress_stream(PyObject *self, PyObject* args) {

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
_decompress_stream(PyObject *self, PyObject* args) {

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


PyMethodDef ReCoDeMethods[] = 
{
	{ "_compress_stream", (PyCFunction)_compress_stream, METH_VARARGS, 0 },
	{ "_decompress_stream", (PyCFunction)_decompress_stream, METH_VARARGS, 0 },
	{0,0,0,0}
};

static struct PyModuleDef PyReCoDe = {
	PyModuleDef_HEAD_INIT,
	"pyReCoDe",   /* name of module */
	"ReCoDe docs", /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
			  or -1 if the module keeps state in global variables. */
	ReCoDeMethods
};

PyMODINIT_FUNC
PyInit_PyReCoDe(void)
{
	return PyModule_Create(&PyReCoDe);
	/*
	PyObject *m;
	if (PyType_Ready(&ReCoDeType) < 0)
		return NULL;

	m = PyModule_Create(&pyrecoder);
	if (m == NULL)
		return NULL;

	Py_INCREF(&ReCoDeType);
	PyModule_AddObject(m, "Recoder", (PyObject *)&ReCoDeType);
	return m;
	*/
}