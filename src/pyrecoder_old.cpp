#include <pyrecoder_include.h>
#include "structmember.h"

PyObject *
makelist(uint8_t *array, size_t size) {
	PyObject *l = PyList_New(size);
	for (size_t i = 0; i != size; ++i) {
		PyList_SET_ITEM(l, i, PyLong_FromUnsignedLong(array[i]));
	}
	return l;
}

//buffer creation function
Py_buffer *
createReadBuffer(int nitems, int itemsize, char *fmt) {
	Py_buffer *py_buf = (Py_buffer *)malloc(sizeof(*py_buf));
	py_buf->obj = NULL;
	py_buf->buf = malloc(nitems * itemsize);
	py_buf->len = nitems * itemsize;
	py_buf->readonly = 0;
	py_buf->itemsize = itemsize;
	py_buf->format = fmt;
	py_buf->ndim = 1;
	py_buf->shape = NULL;
	py_buf->strides = NULL;
	py_buf->suboffsets = NULL;
	py_buf->internal = NULL;
	return py_buf;
}

typedef struct {
	PyObject_HEAD
	FILE *file;
	uint64_t current_frame_index;
	int recode_header;
} RecoderObject;

static void 
ReCoDe_dealloc(RecoderObject *self)
{
	if (self->file != NULL) {
		Py_XDECREF(self->file);
	}
	Py_XDECREF(self->current_frame_index);
	Py_XDECREF(self->recode_header);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
ReCoDe_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	RecoderObject *self;
	self = (RecoderObject *)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->file = NULL;
		self->current_frame_index = 0;
		self->recode_header = 0;
	}
	return (PyObject *)self;
}

static PyMemberDef ReCoDe_members[] = {
	{ "file", T_OBJECT_EX, offsetof(RecoderObject, file), 0, "pointer to file" },
	{ "current_frame_index", T_INT, offsetof(RecoderObject, current_frame_index), 0, "index of next frame to be read" },
	{ "recode_header", T_INT, offsetof(RecoderObject, recode_header), 0, "recode header object" },
	{ NULL }  /* Sentinel */
};

PyObject * 
open_file(RecoderObject* self, PyObject* args) {

	char * filename;
	if (!PyArg_ParseTuple(args, "s", &filename)) {
		return Py_BuildValue("s", "Unable to parse filename argument.");
	}
	self->file = fopen(filename, "r");
	return Py_BuildValue("i", 1);

	//fclose(self->file);
	//return PyLong_FromLongLong(input_value + 1);
	//return makelist(buffer, 10);
	//return Py_BuildValue("b", buffer[0]);
	//return Py_BuildValue("O", &pFile);
}

PyObject *
close_file(RecoderObject* self) {
	fclose(self->file);
	self->file = NULL;
	return Py_BuildValue("i", 1);
}

PyObject *
get_recode_header(RecoderObject* self, PyObject* args) {
	if (self->file == NULL) {
		printf("%s", "File not opened. Use open_file() to open a ReCoDe file before calling get_XX().");
		return Py_BuildValue("i", 0);
	}

	Py_buffer view;
	if (!PyArg_ParseTuple(args, "y*", &view)) {
		printf("%s", "Unable to parse argument: memoryview");
		return Py_BuildValue("i", 0);
	}

	size_t items_read = fread((&view)->buf, 1, 1024, self->file);
	return Py_BuildValue("i", 1);
}

/*
Reads the next frame
*/
PyObject *
get_next_frame(RecoderObject* self) {

	if (self->file == NULL) {
		return Py_BuildValue("s", "File not opened. Use open_file() to open a ReCoDe file before calling get_next_frame().");
	}

	Py_buffer *view = createReadBuffer(10, 1, "B");
	PyObject *mv = PyMemoryView_FromBuffer(view);
	size_t items_read = fread(view->buf, 1, 10, self->file);
	return mv;

}

/*
get_frame(frame_index)
Reads frame frame_index
*/
PyObject * 
get_frame(RecoderObject* self, PyObject* args) {

	Py_buffer view;
	uint64_t frame_index;
	if (!PyArg_ParseTuple(args, "ky*", &frame_index, &view)) {
		return Py_BuildValue("s", "Unable to parse argument frame_index");
	}
	//PyObject *mv = PyMemoryView_FromBuffer(view);
	size_t items_read = fread((&view)->buf, 1, 10, self->file);
	return Py_BuildValue("i",1);

}

/*
get_frames(start_frame_index, num_frames)
reads num_frames frames starting from start_frame_index
*/
PyObject *
get_frames(RecoderObject* self, PyObject* args) {

	uint64_t start_frame_index, num_frames;
	if (!PyArg_ParseTuple(args, "kk", &start_frame_index, &num_frames)) {
		return Py_BuildValue("s", "Unable to parse argument. Expected start_frame_index, num_frames");
	}

	Py_buffer *view = createReadBuffer(10, 1, "B");
	PyObject *mv = PyMemoryView_FromBuffer(view);
	size_t items_read = fread(view->buf, 1, 10, self->file);
	return mv;

}

PyMethodDef ReCoDeMethods[] = 
{
	{ "open_file", (PyCFunction)open_file, METH_VARARGS, 0 },
	{ "close_file", (PyCFunction)close_file, METH_NOARGS, 0 },
	{ "get_recode_header", (PyCFunction)get_recode_header, METH_VARARGS, 0 },
	{ "get_next_frame", (PyCFunction)get_next_frame, METH_VARARGS, 0 },
	{ "get_frame", (PyCFunction)get_frame, METH_VARARGS, 0 },
	{ "get_frames", (PyCFunction)get_frames, METH_VARARGS, 0 },
	{0,0,0,0}
};

static PyTypeObject ReCoDeType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"pyrecoder.Recoder",                            /*tp_name*/
	sizeof(RecoderObject),                          /*tp_basicsize*/
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
	"Recoder object",                               /*tp_doc*/
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


static struct PyModuleDef pyrecoder = {
	PyModuleDef_HEAD_INIT,
	"pyReCoDe",   /* name of module */
	"ReCoDe docs", /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
			  or -1 if the module keeps state in global variables. */
	ReCoDeMethods
};

PyMODINIT_FUNC
PyInit_pyrecoder(void)
{
	//return PyModule_Create(&pyrecoder);
	PyObject *m;
	if (PyType_Ready(&ReCoDeType) < 0)
		return NULL;

	m = PyModule_Create(&pyrecoder);
	if (m == NULL)
		return NULL;

	Py_INCREF(&ReCoDeType);
	PyModule_AddObject(m, "Recoder", (PyObject *)&ReCoDeType);
	return m;
}