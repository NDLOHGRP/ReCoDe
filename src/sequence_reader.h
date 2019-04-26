// Standard (non-extended) header is the default
#define SEQ_HEADER_LEN              1024 

enum eHMetadataFormat {
	H_METADATA_UNKNOWN = 0,
	H_METADATA_BOOL, //bool (1 byte)
	H_METADATA_BYTE, //byte (1 byte)
	H_METADATA_SHORT, //short (2 bytes)
	H_METADATA_USHORT, //unsigned short (2 bytes)
	H_METADATA_INT, //int (4 bytes)
	H_METADATA_UINT, //unsigned int (4 bytes)
	H_METADATA_DOUBLE, //double (8 bytes)
	H_METADATA_STRING, //wchar_t[] (variable)
	H_METADATA_BINARY, //BYTE[] (variable)
	H_METADATA_INT64, //__int64 (8 bytes)
	H_METADATA_UINT64 //unsigned __int64 (8 bytes)
};

enum eHImageFormat {
	H_IMAGE_UNKNOWN = 0, //Unknown format
	H_IMAGE_MONO = 100, //Monochrome Image (LSB)
	H_IMAGE_MONO_BAYER = 101, //Raw Bayer Image (treated as H_IMAGE_MONO)
	H_IMAGE_BGR = 200, //BGR Color Image
	H_IMAGE_PLANAR = 300, //Planar Color Image
	H_IMAGE_RGB = 400, //RGB Color Image
	H_IMAGE_BGRx = 500, //BGRx Color Image
	H_IMAGE_YUV422 = 600, //YUV422
	H_IMAGE_YUV422_20 = 610,
	H_IMAGE_UVY422 = 700, //UVY422
	H_IMAGE_UVY411 = 800, //UVY411
	H_IMAGE_UVY444 = 900, //UVY444
	H_IMAGE_BGR555_PACKED = 905, //PhynxRGB
	H_IMAGE_BGR565_PACKED = 906,
	H_IMAGE_MONO_MSB = 112, //Only for > 8 bit per pixel, MSB align litle endian 10 bit: JIHGFEDC BA000000
	H_IMAGE_MONO_BAYER_MSB = 113, //Only for > 8 bit per pixel, MSB align
	H_IMAGE_MONO_MSB_SWAP = 114, //Only for > 8 bit per pixel, MSB align big endian 10 bit: BA000000 JIHGFEDC
	H_IMAGE_MONO_BAYER_MSB_SWAP = 115, //Only for > 8 bit per pixel, MSB align
	H_IMAGE_BGR10_PPACKED = 123, //Only for 10 bit per pixel, LSB align
	H_IMAGE_BGR10_PPACKED_PHOENIX = 124, //Only for 10 bit per pixel, LSB align, RRRRRRRR RR00GGGG GGGGGGBB BBBBBBBB
	H_IMAGE_RGB10_PPACKED_PHOENIX = 125, //Only for 10 bit per pixel, LSB align, BBBBBBBB BB00GGGG GGGGGGRR RRRRRRRR
	H_IMAGE_MONO_PPACKED = 131, //Only for > 8 bit per pixel, MSB align
	H_IMAGE_MONO_BAYER_PPACKED = 132, //Only for > 8 bit per pixel, MSB align
	H_IMAGE_MONO_PPACKED_8448 = 133, //Only for > 8 bit per pixel, MSB align
	H_IMAGE_MONO_BAYER_PPACKED_8448 = 134,//Only for > 8 bit per pixel, MSB align
	H_IMAGE_GVSP_BGR10V1_PACKED = 135, //Only for 10 bit per pixel(From Gige Vision) BBBBBBBB GGGGGGGG RRRRRRRR 00BBGGRR
	H_IMAGE_GVSP_BGR10V2_PACKED = 136, //Only for 10 bit per pixel(From Gige Vision)00BBBBBB BBGGGGGG GGGGRRRR RRRRRRRR
	H_IMAGE_BASLER_VENDOR_SPECIFIC = 1000,
	H_IMAGE_EURESYS_JPEG = 1001,
	H_IMAGE_ISG_JPEG = 1002
};

enum eHCompression {
	H_COMPRESSION_NONE = 0, //Uncompressed sequence
	H_COMPRESSION_JPEG = 1, //JPEG compressed images
	H_COMPRESSION_RLE = 2, //RLE compressed images (stp3 algo)
	H_COMPRESSION_HUFFMAN = 3, //HUFFMAN compressed images (stp3 algo)
	H_COMPRESSION_LZ = 4, //LZ compressed images (stp3 algo)
	H_COMPRESSION_RLE_FAST = 5, //RLE Fast compressed images (ippwrapper)
	H_COMPRESSION_HUFFMAN_FAST = 6, //HUFFMAN Fast compressed images (ippwrapper)
	H_COMPRESSION_LZ_FAST = 7, //LZ Fast compressed images (ippwrapper)
	H_COMPRESSION_H264 = 8, //H264 compressed images
	H_COMPRESSION_WAVELET = 9, //Wavelet compressed images
	H_COMPRESSION_H264_RGB = 16, //H264 compressed images, RGB after uncompressed (specific from some cameras)
	H_COMPRESSION_JPEG_LOSSESS = 17, //JPeg lossless jpeg compressed SEQ
	H_COMPRESSION_CAYER = 18, //Createc cayer compression.
};

typedef struct {
	unsigned long ImageWidth;			//Image width in pixel
	unsigned long ImageHeight;			//Image height in pixel
	unsigned long ImageBitDepth;		//Image depth in bits (8,16,24,32)
	unsigned long ImageBitDepthReal;	//Precise Image depth (x bits)
	unsigned long ImageSizeBytes;		//Size used to store one image.
	eHImageFormat ImageFormat;			//See formats below
} CImageInfo;

typedef struct {
	long MagicNumber;                   //Always 0xFEED
	uint16_t Name[12];                  //Always "Norpix seq\n"
	long Version;                       //Sequence Header Version
	long HeaderSize;                    //Always 1024
	wchar_t Description[512];           //User description
	CImageInfo ImageInfo;               //Image information
	uint32_t AllocatedFrames;			//Number of frames allocated in the sequence
	uint32_t Origin;					//Should be 0 if not Pre/Post recorded
	uint32_t TrueImageSize;				//Number of bytes between the first pixel of each successive images
	double FrameRate;                   //Suggested Frame rate for playback (in fps)
	long DescriptionFormat;             //The content of "Description" 0-unicode 1-ascii 2-data
	uint64_t ReferenceFrame;            //ReferenceFrame index (0 if none)
	uint64_t FixedSize;                 //Fixed size for compressed sequences (0 if none)
	uint64_t Flags;                     //NorPix reserved flags
	long BayerPattern;                  //Bayer pattern used
	long TimeOffsetUS;                  //Time offset applied to each image timestamp
	long ExtendedHeaderSize;            //Size of the extended header (if any)
	eHCompression CompressionFormat;    //The compression used
	long ReferenceTime;                 //Custom Reference Time (time_t format)
	uint16_t ReferenceTimeMS;           //Custom Reference Time (milliseconds part)
	uint16_t ReferenceTimeUS;           //Custom Reference Time (microseconds part)
	uint64_t H264GOP;                   //Group of Picture value if H264 compression is used
	uint64_t H264Bitrate;               //Bitrate if H264 compression is used
	uint32_t JPEGQualityInfo;			//JPEG Format Quality and lossless
	eHImageFormat H264DecodeFormat;     //H264 decode format
	long long IndexOffset;              //Offset of compression index data
	uint32_t OldestFrameIndex;			//Index of the oldest frame (will be > 0 if loop recording)
	uint32_t BytesAlignment;			//Image alignement in bytes (for uncompressed sequences)
	uint8_t Padding[360];               //Unused bytes, reserved for future uses.
} SEQHeader;


void _parseSEQHeader(uint8_t *headerBytes, SEQHeader* header) {

	memcpy(&header->MagicNumber, &headerBytes[0], 4);
	memcpy(&header->Name, &headerBytes[4], 24);
	memcpy(&header->Version, &headerBytes[28], 4);
	memcpy(&header->HeaderSize, &headerBytes[32], 4);
	memcpy(&header->Description, &headerBytes[36], 512);
	memcpy(&header->ImageInfo.ImageWidth, &headerBytes[548], 4);
	memcpy(&header->ImageInfo.ImageHeight, &headerBytes[552], sizeof(header->ImageInfo.ImageHeight));
	memcpy(&header->ImageInfo.ImageBitDepth, &headerBytes[556], sizeof(header->ImageInfo.ImageBitDepth));
	memcpy(&header->ImageInfo.ImageBitDepthReal, &headerBytes[560], sizeof(header->ImageInfo.ImageBitDepthReal));
	memcpy(&header->ImageInfo.ImageSizeBytes, &headerBytes[564], sizeof(header->ImageInfo.ImageSizeBytes));
	memcpy(&header->ImageInfo.ImageFormat, &headerBytes[568], sizeof(header->ImageInfo.ImageFormat));
	memcpy(&header->AllocatedFrames, &headerBytes[572], sizeof(header->AllocatedFrames));
	memcpy(&header->Origin, &headerBytes[576], sizeof(header->Origin));
	memcpy(&header->TrueImageSize, &headerBytes[580], sizeof(header->TrueImageSize));
	memcpy(&header->FrameRate, &headerBytes[584], sizeof(header->FrameRate));
	memcpy(&header->DescriptionFormat, &headerBytes[592], sizeof(header->DescriptionFormat));
	memcpy(&header->ReferenceFrame, &headerBytes[596], sizeof(header->ReferenceFrame));
	memcpy(&header->FixedSize, &headerBytes[600], sizeof(header->FixedSize));
	memcpy(&header->Flags, &headerBytes[604], sizeof(header->Flags));
	memcpy(&header->BayerPattern, &headerBytes[608], sizeof(header->BayerPattern));
	memcpy(&header->TimeOffsetUS, &headerBytes[612], sizeof(header->TimeOffsetUS));
	memcpy(&header->ExtendedHeaderSize, &headerBytes[616], sizeof(header->ExtendedHeaderSize));
	memcpy(&header->CompressionFormat, &headerBytes[620], sizeof(header->CompressionFormat));
	memcpy(&header->ReferenceTime, &headerBytes[624], sizeof(header->ReferenceTime));
	memcpy(&header->ReferenceTimeMS, &headerBytes[628], sizeof(header->ReferenceTimeMS));
	memcpy(&header->ReferenceTimeUS, &headerBytes[630], sizeof(header->ReferenceTimeUS));
	memcpy(&header->H264GOP, &headerBytes[632], sizeof(header->H264GOP));
	memcpy(&header->H264Bitrate, &headerBytes[636], sizeof(header->H264Bitrate));
	memcpy(&header->JPEGQualityInfo, &headerBytes[640], sizeof(header->JPEGQualityInfo));
	memcpy(&header->H264DecodeFormat, &headerBytes[644], sizeof(header->H264DecodeFormat));
	memcpy(&header->IndexOffset, &headerBytes[648], sizeof(header->IndexOffset));
	memcpy(&header->OldestFrameIndex, &headerBytes[656], sizeof(header->OldestFrameIndex));
	memcpy(&header->BytesAlignment, &headerBytes[660], sizeof(header->BytesAlignment));

}


void getSEQHeader(FILE* fp, SEQHeader** seq_header) {

	uint64_t sz = 1024;
	uint8_t *buffer = (uint8_t *)malloc(sz);

	uint64_t read_count = fread(buffer, sizeof(uint8_t), sz, fp);
	if (read_count != sz) {
		printf("An unexpected error occurred.\nActual and expected read counts do not match. Exiting.");
		exit(0);
	}

	_parseSEQHeader(buffer, *seq_header);
	if ((*seq_header)->MagicNumber != 65261) {
		printf("Invalid Sequence file.");
		exit(1);
	}

}

int _get_frame (FILE* seqFile, uint64_t frame_no, uint16_t *buffer, SEQHeader *header) {
	uint64_t offset = 8192 + (frame_no * header->TrueImageSize);
	uint64_t n_PixelsInFrame = header->ImageInfo.ImageHeight * header->ImageInfo.ImageWidth;
	fseek(seqFile, offset, 0);
	uint64_t read_count = fread(buffer, 2, n_PixelsInFrame, seqFile);
	if (read_count < n_PixelsInFrame) {
		return 0;
	}
	else {
		return 1;
	}
}

void _write_frame(FILE* seqFile, uint64_t frame_no, uint16_t *buffer, SEQHeader *header) {
	uint64_t n_PixelsInFrame = header->ImageInfo.ImageHeight * header->ImageInfo.ImageWidth;
	uint64_t offset = 8192 + (frame_no * header->TrueImageSize);
	fseek(seqFile, offset, 0);
	fwrite(buffer, 2, n_PixelsInFrame, seqFile);
}

uint64_t getSEQFrames(uint8_t process_id, FILE* seqFile, uint16_t *buffer, uint64_t frame_start, uint64_t n_Frames, SEQHeader *header) {

	// Important: this method assumes the data has 16-bit depth
	uint64_t offset;
	uint64_t frame_sz_in_bytes = header->ImageInfo.ImageWidth * header->ImageInfo.ImageHeight * (header->ImageInfo.ImageBitDepth / 8);
	uint64_t n_PixelsInFrame = header->ImageInfo.ImageHeight * header->ImageInfo.ImageWidth;

	// load data
	uint64_t frame_count = 0;

	fseek(seqFile, offset, 0);
	for (uint64_t i = frame_start; i < frame_start + n_Frames; i++) {
		
		offset = 8192 + (i * header->TrueImageSize);
		fseek(seqFile, offset, 0);

		uint64_t read_count = fread(&buffer[frame_count*n_PixelsInFrame], 2, n_PixelsInFrame, seqFile);

		recode_print("RCT %d:Read Count: %d; Num. Frames to Process: %d\n", process_id, read_count, n_Frames);
		if (read_count < n_PixelsInFrame) {
			return frame_count;
		}
		
		frame_count++;
	}
	return frame_count;
}

void setSEQFileData (FILE* seqFile, uint16_t *buffer, uint64_t n_Frames, SEQHeader *header) {
	for (uint64_t i = 0; i < n_Frames; i++) {
		_write_frame(seqFile, i, buffer, header);
	}
}

/*
void loadSEQData(const char* filename, uint16_t *buffer, uint32_t frameStart, uint32_t n_Frames) {
	SEQHeader *seq_header = (SEQHeader *)malloc(sizeof(SEQHeader));
	parseSEQHeader(filename, &seq_header);
	// load data
	FILE *fp = fopen(filename, "rb");
	getSEQFrames(fp, buffer, frameStart, n_Frames, seq_header);
}
*/

void writeSEQHeader(FILE* seqFile, SEQHeader* header) {

	fwrite(&header->MagicNumber, 4, 1, seqFile);
	fwrite(&header->Name, 24, 1, seqFile);
	fwrite(&header->Version, 4, 1, seqFile);
	fwrite(&header->HeaderSize, 4, 1, seqFile);
	fwrite(&header->Description, 512, 1, seqFile);
	fwrite(&header->ImageInfo.ImageWidth, 4, 1, seqFile);
	fwrite(&header->ImageInfo.ImageHeight, 4, 1, seqFile);
	fwrite(&header->ImageInfo.ImageBitDepth, 4, 1, seqFile);
	fwrite(&header->ImageInfo.ImageBitDepthReal, 4, 1, seqFile);
	fwrite(&header->ImageInfo.ImageSizeBytes, 4, 1, seqFile);
	fwrite(&header->ImageInfo.ImageFormat, 4, 1, seqFile);
	fwrite(&header->AllocatedFrames, 4, 1, seqFile);
	fwrite(&header->Origin, 4, 1, seqFile);
	fwrite(&header->TrueImageSize, 4, 1, seqFile);
	fwrite(&header->FrameRate, 8, 1, seqFile);
	fwrite(&header->DescriptionFormat, 4, 1, seqFile);
	fwrite(&header->ReferenceFrame, 4, 1, seqFile);
	fwrite(&header->FixedSize, 4, 1, seqFile);
	fwrite(&header->Flags, 4, 1, seqFile);
	fwrite(&header->BayerPattern, 4, 1, seqFile);
	fwrite(&header->TimeOffsetUS, 4, 1, seqFile);
	fwrite(&header->ExtendedHeaderSize, 4, 1, seqFile);
	fwrite(&header->ReferenceTime, 4, 1, seqFile);
	fwrite(&header->ReferenceTimeMS, 2, 1, seqFile);
	fwrite(&header->ReferenceTimeUS, 2, 1, seqFile);
	fwrite(&header->H264GOP, 4, 1, seqFile);
	fwrite(&header->ExtendedHeaderSize, 4, 1, seqFile);
	fwrite(&header->CompressionFormat, 4, 1, seqFile);
	fwrite(&header->ReferenceTime, 4, 1, seqFile);
	fwrite(&header->ReferenceTimeMS, 4, 1, seqFile);
	fwrite(&header->ReferenceTimeUS, 4, 1, seqFile);
	fwrite(&header->H264GOP, 4, 1, seqFile);
	fwrite(&header->H264Bitrate, 4, 1, seqFile);
	fwrite(&header->JPEGQualityInfo, 4, 1, seqFile);
	fwrite(&header->H264DecodeFormat, 4, 1, seqFile);
	fwrite(&header->IndexOffset, 4, 1, seqFile);
	fwrite(&header->OldestFrameIndex, 4, 1, seqFile);
	fwrite(&header->BytesAlignment, 4, 1, seqFile);

}

void createSEQFiles (uint16_t *data, SEQHeader *seq_header, uint64_t n_Frames, const char* directory, int nFiles) {

	int f;
	for (f = 0; f < nFiles; f++) {
		char part_num[4];
		sprintf(part_num, "%03d", f);
		char* filename = concat(concat(directory, part_num), ".seq");
		FILE *fp = fopen(filename, "wb");
		seq_header->AllocatedFrames = n_Frames - f * 2;
		writeSEQHeader(fp, seq_header);
		setSEQFileData(fp, data, n_Frames - f*2, seq_header);
		fclose(fp);
	}

}
