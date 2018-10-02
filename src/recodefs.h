#define MAX_FILENAME_LENGTH			100

#define RC_MODE_REDUCE_ONLY			0
#define RC_MODE_REDUCE_COMPRESS		1

#define COMPRESSION_LEVEL_FASTEST	0
#define COMPRESSION_LEVEL_BEST		1

#ifdef _WIN32
	#define PATH_SEPARATOR '\\'
#else 
	#define PATH_SEPARATOR '/'
#endif

#define VERBOSITY 0

#ifdef _WIN32
	#define recode_print(format,...) if (VERBOSITY == 1) { printf(format, __VA_ARGS__); }
#else
	#define recode_print(format,args...) if (VERBOSITY == 1) { printf(format, ## args); }
#endif