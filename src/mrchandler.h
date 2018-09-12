// Standard (non-extended) header is the default
#define MRC_HEADER_LEN              1024 

// Data Types for MRC -- IMOD standard
#define MRC_INT8                    0
#define MRC_INT16                   1
#define MRC_FLOAT32                 2
#define MRC_COMPLEX64               4
#define MRC_UINT16                  6

typedef struct {
	
    // Meta-information not included in the standard header
    char*	metaname;
    // Pixelsize isn't directly stored in the header but we do store it for user convience
    float	pixelsize[3]; 
    
    // blosc information
    int32_t	blosc_compressor;
    int		blosc_threads;      // default of zero lets blosc guess the number of threads
    uint8_t	blosc_filter;    	// defaults to 2 for BITSHUFFLE
    size_t	blosc_blocksize;  	// default of zero lets blosc guess the number of threads
    uint8_t blosc_clevel;
    
    // MRC fields
    int32_t	mrcType;
    int32_t	dimensions[3];
    int32_t	nStart[3];
    int32_t	mGrid[3];
    float	cellLen[3];
    float	cellAngle[3];
    int32_t	mapColRowSlice[3];
    
    float	max;
    float	min;
    float	mean;
    
    int32_t	spaceGroup;
    int32_t	extendedHeaderSize;
    
    // MRC2000 fields
    float	origin;
    uint8_t	endian[2];   	// TODO: not implemented yet.
    float	std;
    
    float	voltage;       	// in keV
    float	C3;            	// aka spherical aberration in CEOS formulation
    float	gain;          	// counts/primary electron
	
} MRCHeader;

FILE* void parseMRCHeader(const char* filename, MRCHeader **mrc_header) {
	
}

void getFrame (File* mrcFile, MRCHeader **mrc_header, int frame_number) {
	
}

void getNextFrame (File* mrcFile, MRCHeader **mrc_header) {
	
}

void getFrames (File* mrcFile, MRCHeader **mrc_header, int frame_start, int frame_end) {
	
}