// SequenceReader.cpp : Defines the entry point for the console application.
//
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <stdlib.h>
#include "sequence_reader.h"
//#include "sequence_reader.cpp"

int main()
{
	const char* filename = "D:/cbis/images/20161207/12-06-23.598.seq";
	printf("Made it here 1\n");
	SEQReader *reader = new SEQReader(filename);

	uint64_t buffer_sz = reader->header.ImageInfo.ImageWidth * reader->header.ImageInfo.ImageHeight * (reader->header.ImageInfo.ImageBitDepth / 8);
	uint16_t *frameBuffer = (uint16_t *)malloc(buffer_sz);
	reader->get_frame(10, frameBuffer);
	for (int i = 0; i < 30; i++) {
		reader->get_next_frame(frameBuffer);
	}

	uint16_t *framesBuffer = (uint16_t *)malloc(buffer_sz * 10);
	reader->get_frames(40, 10, framesBuffer);

	delete reader;
	return 0;
}

