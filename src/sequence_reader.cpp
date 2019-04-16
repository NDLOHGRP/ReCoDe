
#include "sequence_reader.hpp"
#include <stdlib.h>

SEQReader::SEQReader() {
    printf("Made it to the constructor");
}

SEQReader::SEQReader(const char *filename) {
    printf("Made it to the constructor.\n");
    fp = fopen (filename, "r");
	if (fp == NULL) {
		printf("File '%s' not found.\n", filename);
	}
	load_header();
}

SEQReader::~SEQReader(void) {
	if (fp != NULL) {
		fclose(fp);
	}
}

void SEQReader::load_header(void) {
	if (fp == NULL) {
		printf("File has to be loaded first before calling load_header().\n");
	} else {
		fread(&header.MagicNumber, 4, 1, fp);
		fread(&header.Name, 2, 12, fp);
		fread(&header.Version, 4, 1, fp);
		fread(&header.HeaderSize, 4, 1, fp);
		//fread(&header.Description, 1, 512, fp);
		fseek(fp, 548, 0);
		fread(&header.ImageInfo.ImageWidth, 4, 1, fp);
		fread(&header.ImageInfo.ImageHeight, 4, 1, fp);
		fread(&header.ImageInfo.ImageBitDepth, 4, 1, fp);
		fread(&header.ImageInfo.ImageBitDepthReal, 4, 1, fp);
		fread(&header.ImageInfo.ImageSizeBytes, 4, 1, fp);
		fread(&header.ImageInfo.ImageFormat, 4, 1, fp);
		fread(&header.AllocatedFrames, 4, 1, fp);
		fread(&header.Origin, 4, 1, fp);
		fread(&header.TrueImageSize, 4, 1, fp);
		fread(&header.FrameRate, 8, 1, fp);
		fread(&header.DescriptionFormat, 4, 1, fp);
		fread(&header.ReferenceFrame, 4, 1, fp);
		fread(&header.FixedSize, 4, 1, fp);
		fread(&header.Flags, 4, 1, fp);
		fread(&header.BayerPattern, 4, 1, fp);
		fread(&header.TimeOffsetUS, 4, 1, fp);
		fread(&header.ExtendedHeaderSize, 4, 1, fp);
		fread(&header.CompressionFormat, 4, 1, fp);
		fread(&header.ReferenceTime, 4, 1, fp);
		fread(&header.ReferenceTimeMS, 2, 1, fp);
		fread(&header.ReferenceTimeU, 2, 1, fp);
		fread(&header.H264GOP, 4, 1, fp);
		fread(&header.H264Bitrate, 4, 1, fp);
		fread(&header.JPEGQualityInfo, 4, 1, fp);
		fread(&header.H264DecodeFormat, 4, 1, fp);
		fread(&header.IndexOffset, 8, 1, fp);
		fread(&header.OldestFrameIndex, 4, 1, fp);
		fread(&header.BytesAlignment, 4, 1, fp);
		fread(&header.Padding, 1, 360, fp);

		this->frame_no = 0;
	}
}

SEQHeader SEQReader::get_header(void) {
	return header;
}

void SEQReader::get_frame(uint64_t frame_no, uint16_t *data) {
	uint64_t offset = 8192 + (frame_no * header.TrueImageSize);
	fseek(fp, offset, 0);
	fread(data, 1, header.ImageInfo.ImageSizeBytes, fp);
	this->frame_no = frame_no;
}

void SEQReader::get_next_frame(uint16_t *data) {
	get_frame(this->frame_no, data);
	this->frame_no++;
}

void SEQReader::get_frames(uint64_t frame_no, uint64_t n_frames, uint16_t *data) {
	uint64_t frame_sz_in_bytes = header.ImageInfo.ImageWidth * header.ImageInfo.ImageHeight * (header.ImageInfo.ImageBitDepth / 8);
	uint64_t pointer = 0;
	this->frame_no = frame_no;
	for (int i = 0; i < n_frames; i++) {
		get_next_frame(&data[pointer]);
		pointer += frame_sz_in_bytes;
	}
}