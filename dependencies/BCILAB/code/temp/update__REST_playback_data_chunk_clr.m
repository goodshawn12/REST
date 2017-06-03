REST_playback_data_range = 1+mod(REST_playback_data.smax:REST_playback_data.smax+size(REST_playback_data_chunk_clr,2)-1,REST_playback_data.buffer_len);
REST_playback_data.marker_pos(:,REST_playback_data_range) = 0;
REST_playback_data.buffer(:,REST_playback_data_range) = REST_playback_data_chunk_clr;
REST_playback_data.smax = REST_playback_data.smax + size(REST_playback_data_chunk_clr,2);