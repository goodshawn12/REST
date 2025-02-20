# REST (Real-time EEG Source-mapping Toolbox): Online ICA and Source Localization
---
REST provides an easy-to-use and interactive platform to perform source-space analysis on live-streaming electroencephalographic (EEG) data in near real-time with visualizations of estimated independent source activations, corresponding source scalp topographies, source locations, as well as other measures.



### To use REST

* To run the code correctly in matlab, the users are supposed to install other ToolBoxs in advance.
  Signal Processing Toolbox 
  Statistics and Machine Learning Toolbox 

* Run the test_script.m to see if REST can be correctly set up in your environment. The script will run REST in playback mode with the demo dataset. If a pipeline designer GUI pops up, just close it so that REST can continue. In the end, you will see the REST visualization GUI pop up that shows the playback data stream, online decomposed independent components, etc. See the Video Tutorials section below for instruction of how to interact with REST interface.

* To work with your own data stream, create a Lab Streaming Layer (LSL) stream for the data you wish REST to work with. To do this, look at the lab streaming layer github at https://github.com/sccn/labstreaminglayer to find how to interface with your data recording setup.

* Open test_script.m to see the example script. Set opts.playback to 0 (false) to enable online streaming mode. If calibration data is defined, the pipeline will be initialized with the data. Otherwise, the pipeline will be initialized with data (default 10-sec) collected on-line. If a headModel is defined, "localize" option in REST will be enable. Type `` help REST`` for more details.

* REST is still in beta version. If you encounter any problem using REST or have any suggestion, please contact us at shh078 at ucsd dot edu.


### Video Tutorials
New Video of REST(2.0) for Real-time Automatic Artifact Rejection! [https://www.youtube.com/watch?v=N3N-rHWexIs].
Previous demo videos: The (first video)[https://www.youtube.com/watch?v=W5mAAS0xxio] is an overview of the layout and basic features of REST. The (second video)[https://www.youtube.com/watch?v=Ud2XNQy1FkE] is an example of how to use REST for real-time eye-component rejection.


### Files description:

- test_script.m		- Example script that shows how to run REST (streaming example data).
- make_headModel.m	- Make your own head model for source localization. Type ``help make_headModel`` for more information.
- /functions/		- All REST functions including online ICA and source localization.
	- /pipeline/		- Contains functions for the processing pipeline
	- /vis/			- Contains functions for visualizations
		- vis_stream_ORICA.m	- Altered version of BCILAB's vis_stream.m for use in REST to show time series.
 - /GUI/			- GUI files and figures.
	- REST.m (.fig)		- Main GUI for visualization and interacting with data stream and sources.
	- selStream.m (.fig)	- selStream GUI allows user to select a LSL stream that REST works with.  
- /data/			- Demo datasets, channel locations, head models, and config files.
- /dependencies/		- Necessary functions from BCILAB, MoBILAB, liblsl and other sources. 

*  flt_orica.m - BCILAB filter function for Online Recursive ICA [2]. Locate in ``\dependencies\BCILAB\code\filters``.


### References
[1] L. Pion-Tonachini, S.-H. Hsu, C.-Y. Chang, T.-P. Jung, S. Makeig, "Automatic Artifact Rejection using the Real-time EEG Source-mapping Toolbox (REST)," in *IEEE EMBS*, 2018.

[2] L. Pion-Tonachini, S.-H. Hsu, S. Makeig, T.-P. Jung, G. Cauwenberghs, "Real-Time EEG Source-Mapping Toolbox (REST): Online ICA and Source Localization," in *IEEE EMBS*, 2015.

[3] S.-H. Hsu, T. R. Mullen, T.-P. Jung and G. Cauwenberghs, "Real-time adaptive EEG source separation using online recursive independent component analysis," *IEEE transactions on neural systems and rehabilitation engineering*, 24(3), pp.309-319, 2016.

2017/07/06 Luca Pion-Tonachini / Sheng-Hsiou Hsu @ SCCN, UCSD
