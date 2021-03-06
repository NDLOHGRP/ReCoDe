from multiprocessing.managers import BaseManager

import numpy as np
from multiprocessing import Process, Manager
import zmq
import time
import json
import threading
from datetime import datetime

class rc:
    REQ_TYPE_QUERY = 0
    REQ_TYPE_COMMAND = 1

    FILE_TYPE_BINARY = 0
    FILE_TYPE_MRC = 1
    FILE_TYPE_SEQ = 2
    FILE_TYPE_OTHER = 255

    STATUS_CODE_BUSY = 0  # Processing a request, can't listen but not dead
    STATUS_CODE_AVAILABLE = 1  # Listening
    STATUS_CODE_ERROR = -1  # Dead due to exception
    STATUS_CODE_NOT_READY = -2  # Hasn't stated yet
    STATUS_CODE_IS_CLOSED = -3  # Is closed
    '''
    Expected status flow: -2, 1, 0, 1, 0, 1, 0, -3
    '''
    STATUS_CODES = [STATUS_CODE_BUSY, STATUS_CODE_AVAILABLE, STATUS_CODE_ERROR, STATUS_CODE_NOT_READY,
                    STATUS_CODE_IS_CLOSED]

    MESSAGE_TYPE_INFO_RESPONSE = 0
    MESSAGE_TYPE_ERROR_RESPONSE = -1
    MESSAGE_TYPE_STATUS_RESPONSE = 1
    MESSAGE_TYPE_ACK_RESPONSE = 2
    MESSAGE_TYPE_REQUEST = 3

    MESSAGE_TYPES = [MESSAGE_TYPE_INFO_RESPONSE, MESSAGE_TYPE_ERROR_RESPONSE, MESSAGE_TYPE_STATUS_RESPONSE,
                     MESSAGE_TYPE_ACK_RESPONSE, MESSAGE_TYPE_REQUEST]

    MESSAGE_TYPE_INFORMAL_DESC = {MESSAGE_TYPE_INFO_RESPONSE: 'INFO', MESSAGE_TYPE_ERROR_RESPONSE: 'INFO',
                                  MESSAGE_TYPE_STATUS_RESPONSE: 'STATUS UPDATE',
                                  MESSAGE_TYPE_ACK_RESPONSE: 'ACKNOWLEDGEMENT',
                                  MESSAGE_TYPE_REQUEST: 'REQUEST'}


class MessageData:

    def __init__(self, session_id, message_type, message='', descriptive_message='', mapped_data=None,
                 source_process_id=None, target_process_id=None):

        self._message_payload = {'session_id': session_id}
        if message_type not in rc.MESSAGE_TYPES:
            raise ValueError('Unknown message type')
        self._message_payload['type'] = message_type
        self._message_payload['message'] = message
        self._message_payload['descriptive_message'] = descriptive_message
        self._message_payload['source_process_id'] = source_process_id
        self._message_payload['target_process_id'] = target_process_id

        if mapped_data is None:
            self._message_payload['mapped_data'] = {}
        else:
            self._message_payload['mapped_data'] = mapped_data

    @property
    def session_id(self):
        return self._message_payload['session_id']

    @property
    def type(self):
        return self._message_payload['type']

    @property
    def message(self):
        return self._message_payload['message']

    @property
    def descriptive_message(self):
        return self._message_payload['descriptive_message']

    @property
    def source_process_id(self):
        return self._message_payload['source_process_id']

    @property
    def target_process_id(self):
        return self._message_payload['target_process_id']

    @property
    def mapped_data(self):
        return self._message_payload['mapped_data']

    def set_mapped_data(self, d):
        self._message_payload['mapped_data'] = d

    def add_mapped_data(self, d):
        self._message_payload['mapped_data'].update(d)

    def from_json(self, json_str):
        self._message_payload = json.loads(json_str)
        return self

    def to_json(self):
        return json.dumps(self._message_payload)

    def print(self):
        print(self._message_payload)


class NodeToken(object):

    def __init__(self, node_id, context, ip_address, server_port, publishing_port):
        self._node_id = node_id
        self._context = context
        self._ip_address = ip_address
        self._server_port = server_port
        self._publishing_port = publishing_port

    @property
    def node_id(self):
        return self._node_id

    @property
    def context(self):
        return self._context

    @property
    def ip_address(self):
        return self._ip_address

    @property
    def server_port(self):
        return self._server_port

    @property
    def publishing_port(self):
        return self._publishing_port


class NodeClient:

    def __init__(self, node_token, node_status=rc.STATUS_CODE_NOT_READY):
        self._node_token = node_token
        self._comm_history = []
        self._server_socket = None
        self._sub_socket = None
        self._is_connected = False

    @property
    def token(self):
        return self._node_token

    def connect(self):
        if self._is_connected:
            return
        try:
            self._server_socket = self._node_token.context.socket(zmq.REQ)
            self._server_socket.connect(self._node_token.ip_address + ":" + str(self._node_token.server_port))
            self._is_connected = True
        except:
            print('Unable to connect to node')
            self._is_connected = False

    def subscribe(self):
        try:
            self._sub_socket = self._node_token.context.socket(zmq.SUB)
            self._sub_socket.connect(self._node_token.ip_address + ":" + str(self._node_token.publishing_port))
        except:
            print('Unable to connect to node')

    def add_to_history(self, message):
        self._comm_history.append(message)

    def send_request(self, md_out):
        self._server_socket.send_json(md_out.to_json())
        s = self._server_socket.recv_json()
        md = MessageData('', 0).from_json(s)
        print('RQM Received:')
        md.print()
        if md.session_id == md_out.session_id:
            if md.mapped_data['req_id'] == md_out.mapped_data['req_id']:
                if md.type == rc.MESSAGE_TYPE_ACK_RESPONSE and md.message == 'ack':
                    return True
        return False

    def close(self):
        if self._is_connected:
            try:
                self._server_socket.close()
                self._is_connected = False
            except Exception as e:
                print(e)


class Logger:

    def __init__(self, token):
        self._node_token = token
        self._socket = None
        self._logs = {}
        self._state = None
        self._session_id = None

    def start(self, session_id, state):
        self._session_id = session_id
        self._state = state
        self._socket = self._node_token.context.socket(zmq.PULL)
        self._socket.bind(self._node_token.ip_address + ":" + str(self._node_token.publishing_port))
        self._state[0] = rc.STATUS_CODE_AVAILABLE
        self.run()

    def run(self):
        self._log(rc.MESSAGE_TYPE_INFO_RESPONSE, 'Logger is up and running')
        while True:
            s = self._socket.recv_json()
            md = MessageData('', 0).from_json(s)
            if md.session_id == self._session_id:
                if md.message == 'close':
                    self._shutdown()
                    self._state[0] = rc.STATUS_CODE_IS_CLOSED
                    self._log(rc.MESSAGE_TYPE_INFO_RESPONSE, 'Logger out')
                    return
                self._push(md)

    def _log(self, severity, text, descriptive_text=''):
        m = MessageData(self._session_id, severity, text,
                        descriptive_message=descriptive_text, source_process_id=self._node_token.node_id,
                        mapped_data={'timestamp': datetime.now().strftime('%m/%d/%Y-%H:%M:%S')})
        self._push(m)

    def _push(self, md):
        if 'timestamp' in md.mapped_data:
            timestamp = md.mapped_data['timestamp']
        else:
            timestamp = datetime.now().strftime("%m/%d/%Y-%H:%M:%S")
        key = timestamp + '_' + str(md.session_id) + '_' + str(md.source_process_id) + '_' + str(md.type)
        self._logs[key] = md
        self._print_msg(timestamp, md)

    @staticmethod
    def _print_msg(timestamp, md):

        if md.source_process_id is not None and md.source_process_id != '':
            if md.source_process_id == -1:
                src = 'Logger'
            elif md.source_process_id == 0:
                src = 'Head Node'
            else:
                src = 'Node ' + str(md.source_process_id)
        else:
            src = ''

        typ = rc.MESSAGE_TYPE_INFORMAL_DESC[md.type]

        print('(' + typ + ')', timestamp, src + ':', md.message, md.descriptive_message,
              '[Session ID:', str(md.session_id) + ']')

    def _shutdown(self):
        self._socket.close()


class ReCoDeWriter:

    def __init__(self):

        self._mode = 'batch'
        self._num_threads = 3
        self._init_params = {}
        self._input_params = {}
        self._calibration_frame = {}

        self._logger_state = None
        self._node_states = None
        self._node_state_update_timestamps = None
        self._context = None
        self._pub_socket = None
        self._node_token = None
        self._session_id = None

        self._max_attempts = 10
        self._max_non_responsive_time = 15

    def run(self):

        if self._mode == 'batch':

            manager = Manager()
            self._logger_state = manager.dict()
            self._node_states = manager.dict()
            self._node_state_update_timestamps = manager.dict()
            self._context = zmq.Context()

            self._session_id = np.random.randint(10000, 11000)
            nodes = []
            node_clients = []

            self._node_token = NodeToken(0, self._context, 'tcp://127.0.0.1', 18534, 28534)

            for i in range(1, self._num_threads+1):
                token = NodeToken(i, self._context, 'tcp://127.0.0.1', 18534 + i, 28534)
                node = NodeClient(token)
                node_clients.append(node)
                self._node_states[token.node_id] = rc.STATUS_CODE_NOT_READY
                self._node_state_update_timestamps[token.node_id] = datetime.now()
                rct = ReCoDeNode(token, self._calibration_frame, self._init_params, self._input_params)
                p = Process(target=rct.run,
                            args=(self._session_id, token,
                                  self._node_states, self._node_state_update_timestamps, self._logger_state))
                nodes.append(p)
                p.start()

            time.sleep(0.1)

            # start logger
            token = NodeToken(-1, self._context, 'tcp://127.0.0.1', -1, 28534)
            logger = Logger(token)
            self._logger_state = {0: rc.STATUS_CODE_NOT_READY}
            logger_process = Process(target=logger.start, args=(self._session_id, self._logger_state))
            logger_process.start()

            time.sleep(0.1)

            # connect to logger
            self._pub_socket = self._context.socket(zmq.PUSH)
            self._pub_socket.connect(self._node_token.ip_address + ":" + str(self._node_token.publishing_port))
            self._log(rc.MESSAGE_TYPE_INFO_RESPONSE, 'Welcome to this session')

            # connect to RC nodes
            for i in range(self._num_threads):
                node_clients[i].connect()

            self._broadcast(node_clients, self._session_id, 1, 'start')
            self._broadcast(node_clients, self._session_id, 2, 'process_file')

            # close RC nodes
            self._broadcast(node_clients, self._session_id, 3, 'close')
            for i in range(self._num_threads):
                node_clients[i].close()
                nodes[i].join()

            # close logger
            self._log(rc.MESSAGE_TYPE_INFO_RESPONSE, 'close')
            logger_process.join()

    def _spawn_replacement_node(self, node_id):
        return

    def _broadcast(self, nodes, session_id, req_id, msg_text):

        _node_acknowledged = np.zeros((len(nodes)), dtype=np.bool)
        _node_attempts = np.zeros((len(nodes)), dtype=np.int)

        n_available_nodes = 0
        for node in nodes:
            if self._node_states[node.token.node_id] != rc.STATUS_CODE_ERROR:
                n_available_nodes += 1

        while np.min(_node_attempts) < self._max_attempts and n_available_nodes > 0:
            for i, node in enumerate(nodes):
                # if rct thread i has not already acknowledged
                if self._node_states[node.token.node_id] != rc.STATUS_CODE_ERROR and not _node_acknowledged[i]:
                    # try (re)sending a message
                    m = MessageData(session_id, rc.MESSAGE_TYPE_REQUEST, msg_text, mapped_data={'req_id': req_id})
                    if self._node_states[node.token.node_id] == rc.STATUS_CODE_AVAILABLE:
                        _node_acknowledged[i] = node.send_request(m)
                        print('RQM sent to server', i, '(port', node.token.server_port, '):')
                        m.print()
                        if not _node_acknowledged[i]:
                            if _node_attempts[i] > self._max_attempts:
                                node.set_status(rc.STATUS_CODE_ERROR)
                            else:
                                _node_attempts[i] += 1
                    else:
                        lct = self._node_state_update_timestamps[node.token.node_id]
                        non_responsive_time = (datetime.now() - lct).total_seconds()
                        if non_responsive_time > self._max_non_responsive_time:
                            print('Head Node: Server', i,
                                  'is non-responsive. (Last contact time:', lct, '). Starting replacement node...')
                            self._spawn_replacement_node(i)
                        else:
                            print('Head Node: Server', i, 'is busy. (Last contact time:', lct, '). Continuing...')

            # check if all nodes have acknowledged
            if np.sum(_node_acknowledged) == len(nodes):
                return

            n_available_nodes = 0
            for node in nodes:
                if self._node_states[node.token.node_id] != rc.STATUS_CODE_ERROR:
                    n_available_nodes += 1

            time.sleep(1)

        return _node_acknowledged, n_available_nodes

    def _log(self, severity, text, descriptive_text=''):
        m = MessageData(self._session_id, severity, text,
                        descriptive_message=descriptive_text, source_process_id=0,
                        mapped_data={'timestamp': datetime.now().strftime('%m/%d/%Y-%H:%M:%S')}).to_json()
        self._pub_socket.send_json(m)


class ReCoDeNode:

    def __init__(self, node_token, calibration_frame, init_params, input_params):
        self._node_token = node_token
        self._pid = node_token.node_id
        self._sport = node_token.server_port
        self._pub_port = node_token.publishing_port
        self._calibration_frame = calibration_frame
        self._init_params = init_params
        self._input_params = input_params

        self._socket = None
        self._pub_socket = None
        self._node_states = None
        self._node_state_update_timestamps = None
        self._logger_state = None
        self._session_id = None

        self._frames = []
        self._timestamps = []

    def _set_state(self, state):
        self._node_states[self._pid] = state
        self._node_state_update_timestamps[self._pid] = datetime.now()

    def _log(self, severity, text, descriptive_text=''):
        m = MessageData(self._session_id, severity, text,
                        descriptive_message=descriptive_text, source_process_id=self._pid,
                        mapped_data={'timestamp': datetime.now().strftime('%m/%d/%Y-%H:%M:%S')}).to_json()
        self._pub_socket.send_json(m)

    def run(self, session_id, token, node_states, node_state_update_timestamps, logger_state):

        self._node_token = token
        self._node_states = node_states
        self._node_state_update_timestamps = node_state_update_timestamps
        self._logger_state = logger_state
        self._session_id = session_id

        # start server
        try:
            context = zmq.Context()

            self._socket = context.socket(zmq.REP)
            self._socket.bind(self._node_token.ip_address + ":" + str(self._sport))

            self._pub_socket = context.socket(zmq.PUSH)
            self._pub_socket.connect(self._node_token.ip_address + ":" + str(self._pub_port))

            self._set_state(rc.STATUS_CODE_AVAILABLE)
            print('Node', self._pid, 'Ready at port', self._sport)

        except Exception as e:
            self._set_state(rc.STATUS_CODE_ERROR)
            print('Node:', self._node_token.node_id, '=>', 'Exception:', e)
            return

        while True:
            print('Node', self._pid, 'Waiting...')
            self._log(rc.MESSAGE_TYPE_INFO_RESPONSE, 'Node ' + str(self._pid) + ' Waiting...')
            s = self._socket.recv_json()
            md = MessageData('', 0).from_json(s)
            print('Node', self._pid, ' Received:')
            md.print()
            if md.session_id == session_id:
                if md.type == rc.MESSAGE_TYPE_REQUEST:

                    if md.message == 'start':
                        self._send_ack(md)
                        self._set_state(rc.STATUS_CODE_BUSY)
                        try:
                            self._start()
                            self._set_state(rc.STATUS_CODE_AVAILABLE)
                        except Exception as e:
                            print(e)
                            self._set_state(rc.STATUS_CODE_ERROR)

                    elif md.message == 'process_file':
                        self._send_ack(md)
                        self._set_state(rc.STATUS_CODE_BUSY)
                        try:
                            self._process_file()
                            self._set_state(rc.STATUS_CODE_AVAILABLE)
                        except Exception as e:
                            print(e)
                            self._set_state(rc.STATUS_CODE_ERROR)

                    elif md.message == 'close':
                        self._send_ack(md)
                        try:
                            self._close(node_states)
                            self._set_state(rc.STATUS_CODE_IS_CLOSED)
                            return
                        except Exception as e:
                            print(e)
                            self._set_state(rc.STATUS_CODE_ERROR)
                            return

                    else:
                        print('Unknown request')
                        self._set_state(rc.STATUS_CODE_ERROR)

    def _send_ack(self, md):
        m = MessageData(
            md.session_id, rc.MESSAGE_TYPE_ACK_RESPONSE, 'ack', mapped_data={'req_id': md.mapped_data['req_id']}
        ).to_json()
        self._socket.send_json(m)

    def _open(self):
        # create output rc part file, serialise recode header
        time.sleep(1)

    def _start(self):
        time.sleep(1)

    def _process_file(self):
        # read source frames
        # reduce-compress and serialize to destination
        # at the end of each run, update q with the frame_ids processed
        self._frames = np.arange(1, 10).tolist()
        self._timestamps = np.arange(1, 10).tolist()
        time.sleep(np.random.randint(5, 10))
        self._log(rc.MESSAGE_TYPE_INFO_RESPONSE, 'Processed 10 frames')
        return

    def _close(self, node_states):
        # update recode header with true frame count, flush remaining data and close file
        self._socket.close()
        return


if __name__ == "__main__":
    writer = ReCoDeWriter()
    writer.run()
