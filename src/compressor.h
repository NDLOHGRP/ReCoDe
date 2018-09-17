

/*
========================================================================================== 
DEFLATE
========================================================================================== 
*/
uint32_t GetMaxCompressedLen_2 (uint32_t  nLenSrc) {
	uint32_t n16kBlocks = (nLenSrc + 16383) / 16384; // round up any fraction of a block
	return (nLenSrc + 6 + (n16kBlocks * 5));
}

float gzip_compress_stream_2 (	uint8_t  *data, 
								uint32_t n_data_bytes, 
								uint8_t  compression_level,
								uint32_t *n_compressed_bytes, 
								uint8_t  *compressedData) {

	//printf("compression_level = %d\n", compression_level);
	
	clock_t p_start = clock();
	
	uint32_t n_dest_bytes = GetMaxCompressedLen_2 (n_data_bytes);

	//unsigned char* pbDst = new unsigned char[nLenDst];  // alloc dest buffer
	//printf("n_data_bytes = %" PRIu32 "\n", n_data_bytes);
	//printf("n_dest_bytes = %" PRIu32 "\n", n_dest_bytes);
	//uint8_t *pbDst = (uint8_t*)calloc(n_dest_bytes, sizeof(uint8_t));
	
	// zlib struct
    z_stream zInfo;
    zInfo.zalloc	= Z_NULL;
    zInfo.zfree 	= Z_NULL;
    zInfo.opaque 	= Z_NULL;
	
    // setup input and compressed output
    zInfo.avail_in 	= n_data_bytes; 			// size of input
	//zInfo.total_in	= n_data_bytes;
    zInfo.avail_out = n_dest_bytes; 			// size of output
	//zInfo.total_out	= n_dest_bytes;
	zInfo.next_in 	= (Bytef *)data; 			// input char array
    zInfo.next_out 	= (Bytef *)compressedData; 	// output char array
	
	// the actual compression work
    int nErr = -1;
	*n_compressed_bytes = 0;
	if (compression_level == 1) {
		nErr = deflateInit(&zInfo, Z_BEST_SPEED); 		// zlib function
	} else if(compression_level == 9) {
		nErr = deflateInit(&zInfo, Z_BEST_COMPRESSION); // zlib function
	}
	if (nErr == Z_OK) {
		nErr = deflate(&zInfo, Z_FINISH);              	// zlib function
		if (nErr == Z_STREAM_END) {
			*n_compressed_bytes = zInfo.total_out;
		}
	}
	deflateEnd(&zInfo);    // zlib function
	
	//*compressedData = pbDst;
	//printf("n_compressed_bytes = %" PRIu32 "\n", *n_compressed_bytes);
	
	//free(pbDst);
	
	//printf("%d\n", *n_compressed_bytes);
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
	
}

float gzip_decompress_stream_2 (uint8_t  compression_level,
								uint8_t  *compressedData, 
								uint8_t  *deCompressedData, 
								uint32_t nCompressedLen, 
								uint32_t nDataLen) {
	
	clock_t p_start = clock();
	
	// zlib struct
    z_stream zInfo;
    zInfo.zalloc	= Z_NULL;
    zInfo.zfree 	= Z_NULL;
    zInfo.opaque 	= Z_NULL;
	
	//zInfo.total_in	= nCompressedLen;
	zInfo.avail_in	= nCompressedLen;
    //zInfo.total_out	= nDataLen;
	zInfo.avail_out	= nDataLen;
    zInfo.next_in	= (Bytef *)compressedData;
    zInfo.next_out	= deCompressedData;

	//printf("%d\n", nCompressedLen);
	//printf("%d\n", nDataLen);
	//printf("%u\n", compressedData[0]);
	
    int nErr, nRet= -1;
    nErr = inflateInit( &zInfo );            // zlib function
    if (nErr == Z_OK) {
		//printf("In 1.\n");
        nErr = inflate(&zInfo, Z_FINISH);   // zlib function
        if (nErr == Z_STREAM_END) {
			//printf("In 2.\n");
            nRet = zInfo.total_out;
        }
    }
    inflateEnd (&zInfo);   					// zlib function
    //printf("De-compressed Size : %d\n", nRet);
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	//printf("gzip_decompress_stream_2 is done.\n");
	return process_time;
}


#if ENABLE_MULTIPLE_COMPRESSIONS

/*
========================================================================================== 
BZIP
========================================================================================== 
*/
float bzip2_compress_stream (	uint8_t  *data, 
								uint32_t n_data_bytes, 
								uint8_t  compression_level,
								uint32_t *n_compressed_bytes, 
								uint8_t  *compressedData) {
									
	clock_t p_start = clock();
	
	uint32_t n_dest_bytes = GetMaxCompressedLen_2 (n_data_bytes);
	
	// zlib struct
    bz_stream bzInfo;
    bzInfo.bzalloc	= Z_NULL;
    bzInfo.bzfree 	= Z_NULL;
    bzInfo.opaque 	= Z_NULL;
	
	bzInfo.avail_in		= n_data_bytes;
    bzInfo.avail_out	= n_dest_bytes;
    
	bzInfo.next_in	= (char *)data;
    bzInfo.next_out	= (char *)compressedData;

	// the actual compression work
    int nErr = -1;
	*n_compressed_bytes = 0;
	
	nErr = BZ2_bzCompressInit(&bzInfo, compression_level, 0, 30); 		// *bz_stream, blockSize100k, verbosity, workFactor
	
	if (nErr == BZ_OK) {
		nErr = BZ2_bzCompress(&bzInfo, BZ_FINISH); 
		if (nErr == BZ_STREAM_END) {
			*n_compressed_bytes = (bzInfo.total_out_hi32 << 32) + bzInfo.total_out_lo32;
		}
	}
	BZ2_bzCompressEnd(&bzInfo);
	
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	printf("n_compressed_bytes = %"PRIu32"\n", *n_compressed_bytes);
	printf("bzip2_compress_stream is done.\n");
	
	return process_time;
	
}


float bzip2_decompress_stream   (uint8_t  compression_level,
								uint8_t  *compressedData, 
								uint8_t  *deCompressedData, 
								uint32_t nCompressedLen, 
								uint32_t nDataLen) {
									
	clock_t p_start = clock();
	
	
    bz_stream bzInfo;
    bzInfo.bzalloc	= Z_NULL;
    bzInfo.bzfree 	= Z_NULL;
    bzInfo.opaque 	= Z_NULL;
	
	bzInfo.avail_in		= nCompressedLen;
	bzInfo.avail_out	= nDataLen;
    bzInfo.next_in		= (Bytef *)compressedData;
    bzInfo.next_out		= deCompressedData;

    int nErr, nRet= -1;
    nErr = BZ2_bzDecompressInit( &bzInfo, 0, 0); 			// *bz_stream, verbosity, small
    if (nErr == BZ_OK) {
		do {
			nErr = BZ2_bzDecompress(&bzInfo);
		} while (nErr != BZ_STREAM_END);
    }
    BZ2_bzDecompressEnd (&bzInfo);
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	return process_time;
}

/*
========================================================================================== 
LZMA
========================================================================================== 
*/
float lzma_compress_stream (	uint8_t  *data, 
								uint32_t n_data_bytes, 
								uint8_t  compression_level,
								uint32_t *n_compressed_bytes, 
								uint8_t  *compressedData) {

								
	clock_t p_start = clock();
	
	uint32_t preset 	= 9;
	lzma_check check 	= LZMA_CHECK_NONE;
	lzma_stream stream 	= LZMA_STREAM_INIT; /* alloc and init lzma_stream struct */
	
	lzma_action action;
	lzma_ret ret_xz;
	
	/* initialize xz encoder */
	ret_xz = lzma_easy_encoder (&stream, preset, check);
	if (ret_xz != LZMA_OK) {
		fprintf (stderr, "lzma_easy_encoder error: %d\n", (int) ret_xz);
		exit(0);
	}

	stream.next_in = data;
	stream.avail_in = n_data_bytes;

	/* out_buf is clean at this point */
	stream.next_out = compressedData;
	stream.avail_out = n_data_bytes;

	/* compress data */
	ret_xz = lzma_code (&stream, LZMA_FINISH);

	if ((ret_xz != LZMA_OK) && (ret_xz != LZMA_STREAM_END)) {
		fprintf (stderr, "lzma_code compress error: %d\n", (int) ret_xz);
		exit(0);
	} else {
		/* write compressed data */
		*n_compressed_bytes = stream.avail_out;	
	}

	lzma_end (&stream);
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	//printf("n_compressed_bytes = %"PRIu32", compresstion_time = %f\n", *n_compressed_bytes, process_time);
	//printf("bzip2_compress_stream is done.\n");
	
	return process_time;
  

}

float lzma_decompress_stream   (uint8_t  compression_level,
								uint8_t  *compressedData, 
								uint8_t  *deCompressedData, 
								uint32_t nCompressedLen, 
								uint32_t nDataLen) {
									
	clock_t p_start = clock();
	
	uint32_t preset 	= 9;
	lzma_check check 	= LZMA_CHECK_NONE;
	lzma_stream stream 	= LZMA_STREAM_INIT; /* alloc and init lzma_stream struct */
	
	lzma_action action;
	lzma_ret ret_xz;
	
	/* initialize xz decoder*/
	ret_xz = lzma_stream_decoder (&stream, UINT64_MAX, LZMA_IGNORE_CHECK);
	if (ret_xz != LZMA_OK) {
		fprintf (stderr, "lzma_stream_decoder error: %d\n", (int) ret_xz);
		exit(0);
	}
	
	stream.next_in 		= compressedData;
	stream.avail_in 	= nCompressedLen;
	stream.next_out 	= deCompressedData;
	stream.avail_out 	= nDataLen;
	
	/* decompress data */
	int count = 0;
	do {
		ret_xz = lzma_code (&stream, LZMA_FINISH);
		fprintf (stderr, "lzma_code: %d, Avail_out: %zu\n", (int) ret_xz, stream.avail_out);
		count++;
	} while (count < 5 && ret_xz == LZMA_OK && ret_xz != LZMA_STREAM_END);
	
	if (ret_xz != LZMA_OK && ret_xz != LZMA_STREAM_END) {
		fprintf (stderr, "lzma_code decompress error: %d\n", (int) ret_xz);
		//exit(0);
	}
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	return process_time;
}

/*
========================================================================================== 
SNAPPY
========================================================================================== 
*/
float snappy_compress_stream (	uint8_t  *data, 
								uint32_t n_data_bytes, 
								uint8_t  compression_level,
								uint32_t *n_compressed_bytes, 
								uint8_t  *compressedData) {

								
	clock_t p_start = clock();


	char* input = (char*) data;
	char* output = (char*) compressedData;
	*n_compressed_bytes = snappy_max_compressed_length(n_data_bytes);
	snappy_status status = snappy_compress (input, n_data_bytes, output, (size_t *)n_compressed_bytes);

	printf("Snappy Compression: Status = %d\n", status);
	printf("Snappy Compression: Compressed Sz = %d\n", *n_compressed_bytes);
	
	if (status != SNAPPY_OK) {
    	fprintf (stderr, "snappy compression error: snappy_status =  %d\n", status);
		exit(0);
  	}

  	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	//printf("n_compressed_bytes = %"PRIu32", compresstion_time = %f\n", *n_compressed_bytes, process_time);
	//printf("bzip2_compress_stream is done.\n");
	
	return process_time;
  

}

float snappy_decompress_stream (uint8_t  compression_level,
								uint8_t  *compressedData, 
								uint8_t  *deCompressedData, 
								uint32_t nCompressedLen, 
								uint32_t nDataLen) {
									
									
	clock_t p_start = clock();

	char* compressed = (char*) compressedData;
	char* uncompressed = (char*) deCompressedData;
	
	size_t expectedDataLength = 0;
	snappy_status status_1 = snappy_uncompressed_length(compressed, nCompressedLen, &expectedDataLength);
	printf("Snappy Status: %d\n", status_1);
	
	status_1 = snappy_validate_compressed_buffer(compressed, (size_t)nCompressedLen);
	printf("Snappy Status: %d\n", status_1);
	
	size_t uncompressedDataLen = 0;
	snappy_status status = snappy_uncompress(compressed, (size_t)nCompressedLen, uncompressed, &uncompressedDataLen);
	
	printf("Snappy: Expected Uncompressed Data Length: %"PRIu32"\n", nDataLen);
	printf("Snappy: Expected Uncompressed Data Length per Snappy: %zu\n", expectedDataLength);
	printf("Snappy: Uncompressed Data Length: %zu\n", &uncompressedDataLen);
	
	if (status != SNAPPY_OK) {
    	fprintf (stderr, "snappy decompression error: snappy_status =  %d\n", status);
		//exit(0);
  	} else if (uncompressedDataLen != (nDataLen)) {
		fprintf (stderr, "Snappy error: Decompressed size does not match expected size.\n");
		exit(0);
	}
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	return process_time;
}

/*
========================================================================================== 
LZ4
========================================================================================== 
*/
float lz4_compress_stream (	uint8_t  *data, 
							uint32_t n_data_bytes, 
							uint8_t  compression_level,
							uint32_t *n_compressed_bytes, 
							uint8_t  *compressedData) {

								
	clock_t p_start = clock();


	const int max_dst_size = LZ4_compressBound(n_data_bytes);
	char* input = (char*) data;
	char* output = (char*) compressedData;
	*n_compressed_bytes = LZ4_compress_default(input, output, n_data_bytes, max_dst_size);
	if (*n_compressed_bytes < 0) {
		fprintf (stderr, "LZ4 failed with exit code %d\n", *n_compressed_bytes);
		exit(0);
	} else if (*n_compressed_bytes == 0) {
		fprintf (stderr, "LZ4 failed: Destination buffer too small.\n");
		exit(0);
	}

  	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	//printf("n_compressed_bytes = %"PRIu32", compresstion_time = %f\n", *n_compressed_bytes, process_time);
	//printf("bzip2_compress_stream is done.\n");
	
	return process_time;
  

}

float lz4_decompress_stream   (	uint8_t  compression_level,
								uint8_t  *compressedData, 
								uint8_t  *deCompressedData, 
								uint32_t nCompressedLen, 
								uint32_t nDataLen) {
									
	clock_t p_start = clock();
	
	char* compressed_data = (char*) compressedData;
	char* regen_buffer = (char*) deCompressedData;
	const int decompressed_size = LZ4_decompress_safe(compressed_data, regen_buffer, nCompressedLen, nDataLen);
	if (decompressed_size < 0) {
		fprintf (stderr, "LZ4 failed with exit code %d\n", decompressed_size);
		exit(0);
	} else if (decompressed_size == 0) {
		fprintf (stderr, "LZ4 failed: Unknown reason.\n");
		exit(0);
	} else if (decompressed_size != nDataLen) {
		fprintf (stderr, "LZ4 failed: Decompressed size does not match expected size.\n");
		exit(0);
	}
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	return process_time;
}

#endif


#if ENABLE_MULTIPLE_COMPRESSIONS

/*
========================================================================================== 
Main Compression Entry Point
========================================================================================== 
*/
float compress_stream (	uint8_t  compression_scheme,
						uint8_t  compression_level,
						uint8_t  *data, 
						uint32_t n_data_bytes, 
						uint32_t *n_compressed_bytes, 
						uint8_t  *compressedData) {
	
	switch (compression_scheme) {

		case 0:
			return gzip_compress_stream_2 (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

		case 1:
			return bzip2_compress_stream (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

		case 2:
			return lzma_compress_stream (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

		case 3:
			return snappy_compress_stream (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

		case 4:
			return lz4_compress_stream (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

		default:
			return gzip_compress_stream_2 (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

	}

}


/*
========================================================================================== 
Main Decompression Entry Point
========================================================================================== 
*/
float decompress_stream (	uint8_t  compression_scheme,
							uint8_t  compression_level,
							uint8_t  *compressedData, 
							uint8_t  *deCompressedData, 
							uint32_t nCompressedLen, 
							uint32_t nDataLen) {

	switch(compression_scheme) {
		
		case 0:
			return gzip_decompress_stream_2 (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;

		case 1:
			return bzip2_decompress_stream (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;

		case 2:
			return lzma_decompress_stream (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;

		case 3:
			return snappy_decompress_stream (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;

		case 4:
			return lz4_decompress_stream (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;

		default:
			return gzip_decompress_stream_2 (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;
	}

}

#else
	
/*
========================================================================================== 
Main Compression Entry Point
========================================================================================== 
*/
float compress_stream (	uint8_t  compression_scheme,
						uint8_t  compression_level,
						uint8_t  *data, 
						uint32_t n_data_bytes, 
						uint32_t *n_compressed_bytes, 
						uint8_t  *compressedData) {
	
	switch (compression_scheme) {

		case 0:
			return gzip_compress_stream_2 (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

		default:
			return gzip_compress_stream_2 (data, n_data_bytes, compression_level, n_compressed_bytes, compressedData);
			break;

	}

}


/*
========================================================================================== 
Main Decompression Entry Point
========================================================================================== 
*/
float decompress_stream (	uint8_t  compression_scheme,
							uint8_t  compression_level,
							uint8_t  *compressedData, 
							uint8_t  *deCompressedData, 
							uint32_t nCompressedLen, 
							uint32_t nDataLen) {

	switch(compression_scheme) {
		
		case 0:
			return gzip_decompress_stream_2 (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;

		default:
			return gzip_decompress_stream_2 (compression_level, compressedData, deCompressedData, nCompressedLen, nDataLen);
			break;
	}

}

#endif

