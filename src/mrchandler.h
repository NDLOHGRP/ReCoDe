// Standard (non-extended) header is the default
#define MRC_HEADER_LEN              1024 

// Data Types for MRC -- IMOD standard
#define MRC_INT8                    0
#define MRC_INT16                   1
#define MRC_FLOAT32                 2
#define MRC_COMPLEX64               4
#define MRC_UINT16                  6

typedef struct {
	
    // MRC fields
    int32_t	dimensions[3];
    int32_t	mrcType;
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

int _parseMRCHeader(uint8_t *headerBytes, MRCHeader* header) {

	memcpy(&header->dimensions, &headerBytes[0], sizeof(header->dimensions));
	memcpy(&header->mrcType, &headerBytes[12], sizeof(header->mrcType));

	/*
	if( header->mrcType >= MRC_COMP_RATIO )
	{   // Check if the mrc type/mode is >= 1000, which indicates a compressed type
	header->blosc_compressor = header->mrcType / MRC_COMP_RATIO;
	header->mrcType = header->mrcType % MRC_COMP_RATIO;
	}
	else
	{
	header->blosc_compressor = BLOSC_COMPRESSOR_NONE;
	}
	*/

	memcpy(&header->nStart, &headerBytes[16], sizeof(header->nStart));
	memcpy(&header->mGrid, &headerBytes[28], sizeof(header->mGrid));
	memcpy(&header->cellLen, &headerBytes[40], sizeof(header->cellLen));
	memcpy(&header->cellAngle, &headerBytes[52], sizeof(header->cellAngle));
	memcpy(&header->mapColRowSlice, &headerBytes[64], sizeof(header->mapColRowSlice));

	memcpy(&header->min, &headerBytes[76], sizeof(header->min));
	memcpy(&header->max, &headerBytes[80], sizeof(header->max));
	memcpy(&header->mean, &headerBytes[84], sizeof(header->mean));

	memcpy(&header->spaceGroup, &headerBytes[88], sizeof(header->spaceGroup));
	memcpy(&header->extendedHeaderSize, &headerBytes[92], sizeof(header->extendedHeaderSize));

	// MRC2000 fields
	memcpy(&header->origin, &headerBytes[196], sizeof(header->origin));
	memcpy(&header->endian, &headerBytes[212], sizeof(header->endian));
	memcpy(&header->std, &headerBytes[216], sizeof(header->std));

	// MRCZ fields
	memcpy(&header->voltage, &headerBytes[132], sizeof(header->voltage));
	memcpy(&header->C3, &headerBytes[136], sizeof(header->C3));
	memcpy(&header->gain, &headerBytes[140], sizeof(header->gain));

	printf("Dimensions: %i %i %i\n", header->dimensions[0], header->dimensions[1], header->dimensions[2]);
	printf("MRC type: %i\n", header->mrcType);
	printf("Cell length: %f, %f, %f\n", header->cellLen[0], header->cellLen[1], header->cellLen[2]);
	printf("Min: %f, Max: %f, Mean: %f\n", header->min, header->max, header->mean);
	printf("Extended header size in bytes: %i\n", header->extendedHeaderSize);

	return 0;
}


FILE* parseMRCHeader(const char* filename, MRCHeader **mrc_header) {
    
    uint64_t sz = 256;
    uint8_t *buffer = (uint8_t *)malloc(sz);
    
    FILE *fp = fopen(filename, "rb");
    uint64_t read_count = fread(buffer, sizeof(uint8_t), sz, fp);
    fclose(fp);
    
	if (read_count != sz) {
		printf("An unexpected error occurred.\nActual and expected read counts do not match. Exiting.");
		exit(0);
	}
    
    _parseMRCHeader(buffer, *mrc_header);
    
	return fp;
}

/*
void getFrame (File* mrcFile, MRCHeader **mrc_header, int frame_number) {
	
}

void getNextFrame (File* mrcFile, MRCHeader **mrc_header) {
	
}
*/
void getMRCFrames (FILE* mrcFile, uint16_t *buffer, uint64_t frame_start, uint64_t n_Frames, DataSize d) {
	
	// Important: this method assumes the data has 16-bit depth
	
	// load data
	size_t seek_count = fseek(mrcFile, MRC_HEADER_LEN + (d.nx*d.ny*frame_start), SEEK_SET);

	uint64_t n_pixels_in_frame = d.nx*d.ny;
	uint64_t chunk_size_frames = 200;
	uint64_t chunk_size_pixels = chunk_size_frames * n_pixels_in_frame;
	uint64_t iter_c = floor((uint64_t)n_Frames / chunk_size_frames);
	uint64_t rem_sz = (uint64_t)n_Frames % chunk_size_frames;
	uint64_t i;
	uint64_t read_count = 0;
	uint64_t n_elem = 0;

	recode_print("nx: %lu\n", d.nx);
	recode_print("ny: %lu\n", d.ny);
	recode_print("nf: %"PRIu64"\n", n_Frames);
	recode_print("iter_c: %"PRIu64"\n", iter_c);

	for (i = 0; i < iter_c; i++) {
		n_elem = i*chunk_size_pixels;
		read_count += fread(buffer + n_elem, sizeof(uint16_t), chunk_size_pixels, mrcFile);
	}
	n_elem = i*chunk_size_pixels;
	read_count += fread(buffer + i*chunk_size_pixels, sizeof(uint16_t), n_pixels_in_frame*rem_sz, mrcFile);
	fclose(mrcFile);

	recode_print("Seek Count: %"PRIu64"\n", seek_count);
	recode_print("Read Count: %"PRIu64"\n", read_count);

	uint64_t expected_read_count = d.nx*d.ny*n_Frames;
	recode_print("Expected Read Count: %"PRIu64"\n", expected_read_count);

	if (read_count != expected_read_count) {
		printf("An unexpected error occurred.\nActual and expected read counts do not match. Exiting.");
		exit(0);
	}

}

int getBitDepth(MRCHeader* header) {
    if( header->mrcType == MRC_INT8 ) {
        return 8;
    } else if( header->mrcType == MRC_INT16 ) {
        return 16;
    } else if( header->mrcType == MRC_FLOAT32 ) {
        return 32;
    } else if( header->mrcType == MRC_COMPLEX64 ) {
        return 64;
    } else if( header->mrcType == MRC_UINT16 ) {
        return 16;
    } else {
		printf("Error in getBitDepth() in mrchandler: unknown MRC type %d.", header->mrcType);
		exit(0);
    }
}
