
#include <inttypes.h>
#include <stdio.h>

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
    H_COMPRESSION_NONE=0, //Uncompressed sequence
    H_COMPRESSION_JPEG=1, //JPEG compressed images
    H_COMPRESSION_RLE=2, //RLE compressed images (stp3 algo)
    H_COMPRESSION_HUFFMAN=3, //HUFFMAN compressed images (stp3 algo)
    H_COMPRESSION_LZ=4, //LZ compressed images (stp3 algo)
    H_COMPRESSION_RLE_FAST=5, //RLE Fast compressed images (ippwrapper)
    H_COMPRESSION_HUFFMAN_FAST=6, //HUFFMAN Fast compressed images (ippwrapper)
    H_COMPRESSION_LZ_FAST=7, //LZ Fast compressed images (ippwrapper)
    H_COMPRESSION_H264=8, //H264 compressed images
    H_COMPRESSION_WAVELET=9, //Wavelet compressed images
    H_COMPRESSION_H264_RGB = 16, //H264 compressed images, RGB after uncompressed (specific from some cameras)
    H_COMPRESSION_JPEG_LOSSESS = 17, //JPeg lossless jpeg compressed SEQ
    H_COMPRESSION_CAYER = 18, //Createc cayer compression.
};

struct CImageInfo {
    unsigned long ImageWidth; //Image width in pixel
	unsigned long ImageHeight; //Image height in pixel
	unsigned long ImageBitDepth; //Image depth in bits (8,16,24,32)
	unsigned long ImageBitDepthReal;//Precise Image depth (x bits)
	unsigned long ImageSizeBytes; //Size used to store one image.
    eHImageFormat ImageFormat; //See formats below
};

struct SEQHeader {
    long MagicNumber;                   //Always 0xFEED
    uint16_t Name[12];                  //Always "Norpix seq\n"
    long Version;                       //Sequence Header Version
    long HeaderSize;                    //Always 1024
	wchar_t Description[512];           //User description
    CImageInfo ImageInfo;               //Image information
	uint32_t AllocatedFrames;      //Number of frames allocated in the sequence
	uint32_t Origin;               //Should be 0 if not Pre/Post recorded
	uint32_t TrueImageSize;        //Number of bytes between the first pixel of each successive images
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
    uint16_t ReferenceTimeU;            //Custom Reference Time (microseconds part)
    uint64_t H264GOP;                   //Group of Picture value if H264 compression is used
    uint64_t H264Bitrate;               //Bitrate if H264 compression is used
	uint32_t JPEGQualityInfo;      //JPEG Format Quality and lossless
    eHImageFormat H264DecodeFormat;     //H264 decode format
    long long IndexOffset;              //Offset of compression index data
	uint32_t OldestFrameIndex;     //Index of the oldest frame (will be > 0 if loop recording)
	uint32_t BytesAlignment;       //Image alignement in bytes (for uncompressed sequences)
    uint8_t Padding[360];               //Unused bytes, reserved for future uses.
};

class SEQReader {
   public:
      SEQHeader header;
      SEQHeader get_header(void);
      void get_next_frame(uint16_t *data);
      void get_frame(uint64_t frame_no, uint16_t *data);
	  void get_frames(uint64_t frame_no, uint64_t n_frames, uint16_t *data);
      void get_frame_time(uint64_t frame_no);
      SEQReader();
      SEQReader(const char *filename); // loads header
      ~SEQReader(void);
   private:
      void load_header( void );
      FILE *fp;
      uint64_t frame_no;
      uint64_t frames_in_buffer;
};