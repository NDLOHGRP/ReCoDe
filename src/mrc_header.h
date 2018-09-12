typedef struct {
	
	/*
	Frame dimensions
	*/
	unsigned short	nx;
	unsigned short	ny;
	
	/*
	Number of frames in the dataset
	*/
	unsigned long	nz;
	
	/*
	The bit depth used to store pixel intensity values. Only used for L1 and L2. Ignored in L3 and L4.
	*/
	unsigned char	bit_depth;
	
	
	/*
	The entire raw data of the MRC header as one chunk
	*/
	unsigned char	data[1024];
	
	
} MRCHeader;

/*
parse_mrc_header (FILE *fp, MRCHeader **mrc_header) {
	fread(&mrc_header->data, sizeof(1024), 1, fp);
}


serialize_mrc_header (FILE *fp, MRCHeader *mrc_header) {
	fwrite (mrc_header->data, sizeof(1024), 1, fp);
}
*/