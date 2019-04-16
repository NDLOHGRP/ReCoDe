#define VERSION_MAJOR				0
#define VERSION_MINOR				1

#define MAX_FILENAME_LENGTH			100

#define RC_MODE_REDUCE_ONLY			0
#define RC_MODE_REDUCE_COMPRESS		1

#define COMPRESSION_LEVEL_FASTEST	0
#define COMPRESSION_LEVEL_BEST		1

#define SOURCE_FILE_TYPE_BINARY		0
#define SOURCE_FILE_TYPE_MRC		1
#define SOURCE_FILE_TYPE_SEQUENCE	2
#define SOURCE_FILE_TYPE_OTHER		3

#define SOURCE_HEADER_POSITION_PRE	0
#define SOURCE_HEADER_POSITION_POST	1

#ifdef _WIN32
	#define PATH_SEPARATOR '\\'
#else 
	#define PATH_SEPARATOR '/'
#endif

int VERBOSITY = 1;

#ifdef _WIN32
	#define recode_print(format,...) if (VERBOSITY == 1) { printf(format, __VA_ARGS__); }
#else
	#define recode_print(format,args...) if (VERBOSITY == 1) { printf(format, ## args); }
#endif

#define	 DTYPE_USHORT	16

#define SetBit(A,k)		( A[(k/8)] |= (1 << (k%8)) )
#define ClearBit(A,k)	( A[(k/8)] &= ~(1 << (k%8)) )
#define CheckBit(A,k)	( A[(k/8)] & (1 << (k%8)) )

typedef struct {
	uint64_t	nx;
	uint64_t	ny;
	uint64_t	nz;
	uint64_t	dtype;
} DataSize;

typedef struct {

	unsigned char	reduction_level;

	/*
	rc_operation_mode can be 0 => reduction only or 1 => reduction and compression,
	there is also a compression only mode, rc_operation_mode = 2 but this is used only for diagnostics
	Default is rc_operation_mode = 1.
	*/
	unsigned char	rc_operation_mode;

	/*
	For the (ij)th pixel the threshold used for dark subtraction is estimated as:
	the maximum value of the (ij)th pixel across the z frames of the dark noise dataset + dark_threshold_epsilon
	*/
	unsigned short	dark_threshold_epsilon;

	/*
	target bit depth after reduction
	if not provided, same as source_bit_depth
	*/
	unsigned char	bit_depth;

	/*
	The bit depth of pixel intensity values in the source file.
	If source file is MRC the bit depth as per MRC header is used
	If source file is not MRC this value must be explicitly provided
	*/
	unsigned char	source_bit_depth;

	/*
	Frame dimensions
	If source file is MRC, these parameters are read from MRC header, else user provided values are used
	*/
	unsigned short	num_cols;
	unsigned short	num_rows;

	/*
	Number of frames to extract from dataset i
	If source file is MRC and the user provided value is -1, the number of frames as per MRC header is used (i.e. all frames are used)
	If source file is MRC and the user provided value is not -1, the user specified number of frames are read
	If source file is not MRC this value has to be explicitly provided
	*/
	unsigned long	num_frames;

	/*
	Number of initial frames to skip in the original dataset
	*/
	unsigned long	frame_offset;

	/*
	Number of frames to be used for dark noise estimation
	If source file is MRC and the user provided value is -1, the number of frames as per MRC header is used (i.e. all frames are used)
	If source file is MRC and the user provided value is not -1, the user specified number of frames are read
	If source file is not MRC this value has to be explicitly provided
	*/
	unsigned long	num_dark_frames;

	/*
	Number of initial frames to skip in the dark noise data set
	*/
	unsigned long	dark_frame_offset;

	unsigned char	keep_part_files;
	unsigned char	num_threads;
	unsigned char	L2_statistics;		// 0 = None (default for L1, L3 and L4), 1 = Max, 2 = Mean
	unsigned char	L4_centroiding;		// 0 = None (default for L1, L2 and L3), 1 = Weighted Centroids, 1 = Max. Pixel, 2 = Unweighted Centroids
	unsigned char	compression_scheme;	// 0 = none, 1 = gzip, 2 = bzip, 3 = lzma (For forward compatibility only. Currently only gzip is supported. The parameter value will be ignored and automatically set to 1.)
	unsigned char	compression_level;	// gzip / bzip compression level (For forward compatibility only. Currently only gzip level 0 is supported. The parameter value will be ignored and automatically set to 0.)
	unsigned char	keep_dark_data;

	unsigned char	source_file_type;
	unsigned short	source_header_length;
	unsigned char	dark_file_type;
	unsigned short	dark_header_length;

} InputParams;


typedef struct {
	char*		run_name;
	char*		start_time;
	uint16_t	thread_id;
	uint64_t	num_writes;
	uint64_t	wait_time;
	uint64_t	reduction_time;
	uint64_t	compression_time;
	uint64_t	reduced_bytes;
	uint64_t	compressed_bytes;
	uint64_t	decompression_time;
	uint64_t	total_time;
} RunStats;
