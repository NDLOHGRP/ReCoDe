import numpy as np
from .params import InitParams, InputParams

class ReCoDeHeader():

    '''
    ToDo:
    Original File Header Position: should be always 0 default'
    Dark Frame Offset and Number of Dark Frames Used for Calibration should be removed
    
    '''

    def __init__(self):
        self._rc_header = {}
        self._rc_header_field_defs = []
        self._rc_header_field_defs.append({'name': 'uid', 'bytes': 8, 'dtype': np.uint64})
        self._rc_header_field_defs.append({'name': 'version_major', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'version_minor', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'reduction_level', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'recode_operation_mode', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'bit_depth', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'nx', 'bytes': 2, 'dtype': np.uint16})
        self._rc_header_field_defs.append({'name': 'ny', 'bytes': 2, 'dtype': np.uint16})
        self._rc_header_field_defs.append({'name': 'nz', 'bytes': 4, 'dtype': np.uint32})
        self._rc_header_field_defs.append({'name': 'L2_statistics', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'L4_centroiding', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'compression_scheme', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'compression_level', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'source_file_type', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'source_header_length', 'bytes': 2, 'dtype': np.uint16})
        self._rc_header_field_defs.append({'name': 'source_header_position', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'source_file_name', 'bytes': 100, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'dark_file_name', 'bytes': 100, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'dark_threshold_epsilon', 'bytes': 2, 'dtype': np.uint16})
        self._rc_header_field_defs.append({'name': 'has_dark_data', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'frame_offset', 'bytes': 4, 'dtype': np.uint32})
        self._rc_header_field_defs.append({'name': 'dark_frame_offset', 'bytes': 4, 'dtype': np.uint32})
        self._rc_header_field_defs.append({'name': 'num_dark_frames', 'bytes': 4, 'dtype': np.uint32})
        self._rc_header_field_defs.append({'name': 'source_bit_depth', 'bytes': 1, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'checksum', 'bytes': 32, 'dtype': np.uint8})
        self._rc_header_field_defs.append({'name': 'futures', 'bytes': 44, 'dtype': np.uint8})

        self._rc_header_length = 321

    def create (self, init_params, input_params):
        self._rc_header['uid'] = 158966344846346
        self._rc_header['version_major'] = 0
        self._rc_header['version_minor'] = 1
        self._rc_header['reduction_level'] = input_params.reduction_level
        self._rc_header['recode_operation_mode'] = input_params.recode_operation_mode
        self._rc_header['bit_depth'] = input_params.bit_depth
        self._rc_header['nx'] = input_params.nx
        self._rc_header['ny'] = input_params.ny
        self._rc_header['nz'] = input_params.nz
        self._rc_header['L2_statistics'] = input_params.L2_statistics
        self._rc_header['L4_centroiding'] = input_params.L4_centroiding
        self._rc_header['compression_scheme'] = input_params.compression_scheme
        self._rc_header['compression_level'] = input_params.compression_level
        self._rc_header['source_file_type'] = input_params.source_file_type
        self._rc_header['source_header_length'] = input_params.source_header_length
        self._rc_header['source_header_position'] = 0
        self._rc_header['source_file_name'] = init_params.image_filename
        self._rc_header['dark_file_name'] = init_params.dark_filename
        self._rc_header['dark_threshold_epsilon'] = input_params.dark_threshold_epsilon
        self._rc_header['has_dark_data'] = input_params.has_dark_data
        self._rc_header['frame_offset'] = input_params.frame_offset
        self._rc_header['dark_frame_offset'] = input_params.dark_frame_offset
        self._rc_header['num_dark_frames'] = input_params.num_dark_frames
        self._rc_header['source_bit_depth'] = input_params.source_bit_depth
        self._rc_header['futures'] = np.array((255), dtype=np.uint8)

    def load(self, rc_filename):
        assert rc_filename != '', 'ReCoDe filename missing'
        with open(rc_filename, 'rb') as fp:
            for field in self._rc_header_field_defs:
                b = fp.read(field['bytes'])
                value = np.frombuffer(b, dtype=field['dtype'])
                if field['name'] is 'dark_file_name':
                    formatted_value = self._to_string(value)
                elif field['name'] is 'source_file_name':
                    formatted_value = self._to_string(value)
                else:
                    if len(value) == 1:
                        formatted_value = value[0]
                    else:
                        formatted_value = value
                self._rc_header[field['name']] = formatted_value

    def serialize(self, rc_filename):
        assert rc_filename != '', 'ReCoDe filename missing'
        with open(rc_filename, 'wb') as fp:
            for field in self._rc_header_field_defs:
                fp.write(self._rc_header[field['name']])

    def skip_header(self, rc_fp):
        rc_fp.seek(self._rc_header_length)
        return rc_fp

    def print(self):
        for field in self._rc_header_field_defs:
            print(field['name'], '=', self._rc_header[field['name']])
            '''
            if field['name'] is 'dark_file_name':
                print(field['name'], '=', self._to_string(self._rc_header[field['name']]))
            elif field['name'] is 'source_file_name':
                print(field['name'], '=', self._to_string(self._rc_header[field['name']]))
            else:
                if len(self._rc_header[field['name']]) == 1:
                    print(field['name'], '=', self._rc_header[field['name']][0])
                else:
                    print(field['name'], '=', self._rc_header[field['name']])
            '''
            
    def validate(self):
        # check that no fields are present in header, does not check the validity of their values
        for field in self._rc_header_field_defs:
            if field['name'] not in self._rc_header:
                return False
        return True

    def _load_metadata(self, rc_filename):
        assert rc_filename != '', 'ReCoDe filename missing'
        # nCompressedSize_BinaryImage, nCompressedSize_Pixvals, bytesRequiredForPacking
        with open(rc_filename, 'wb') as fp:
            buffer = fp.read(self._rc_header['nz']*3*4)
            values = np.frombuffer(buffer, dtype=np.uint32)
            self._metadata = np.reshape(values, [self._rc_header['nz'],3])
        print(self._metadata[:100,:])

    def _to_string(self, arr):
        return ''.join([chr(x) for x in arr])


if __name__== "__main__":

    rc_header = ReCoDeHeader()
    rc_header.load('D:/cbis/GitHub/ReCoDe/scratch/400fps_dose_43.rc1')
    rc_header.print()