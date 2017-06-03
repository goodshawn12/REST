REST_test_data_range = 1+mod(REST_test_data.smax:REST_test_data.smax+size(REST_test_data_chunk_clr,2)-1,REST_test_data.buffer_len);
REST_test_data.marker_pos(:,REST_test_data_range) = 0;
REST_test_data.buffer(:,REST_test_data_range) = REST_test_data_chunk_clr;
REST_test_data.smax = REST_test_data.smax + size(REST_test_data_chunk_clr,2);