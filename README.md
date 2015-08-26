# REST (Real-time EEG Source-mapping Toolbox): 
Online ICA and Source Localization
---
REST provides an easy-to-use and interactive platform to perform source-space analysis on live-streaming electroencephalographic (EEG) data in near real-time with visualizations of estimated independent source activations, corresponding source scalp topographies, source locations, as well as other measures.



### To use REST

* Make sure the BCILAB is called before running the main script. The latest versions of BCILAB can be found at github.com/sccn/BCILAB. When downloaded, go to bcilab directory and type ``bcilab`` to initiate it. Alternatively, you can type out the path to the BCILAB directory into the  txt file *bcilab_path.txt* and REST will automatically start bcilab if necessary (recommended).

* Create a Lab Streaming Layer (LSL) stream for the data you wish REST to work with. To do this, look at the lab streaming layer wiki at https://code.google.com/p/labstreaminglayer/w/list to find how to interface with your data recording setup. LSL is included in BCILAB but can also be found at https://code.google.com/p/labstreaminglayer/. 

* Start REST with the desired settings. Type `` help REST`` to see options.

### Files description:

- test_script.m		- Example script that shows how to run REST (streaming example data).
- make_headModel.m	- Make your own head model for source localization. Type ``help make_headModel`` for more information.
- bcilab_path.txt	- Define path to BCILAB.
- /functions/		- All REST functions including online ICA and source localization.
	- /pipeline/		- Contains functions for the processing pipeline
		- flt_orica.m		- BCILAB filter function for ORICA. Calls dynamicWhitening.m and dynamicOrica.m
	- /vis/			- Contains functions for visualizations
		- vis_stream_ORICA.m	- Altered version of BCILAB's vis_stream.m for use in REST to show time series.
 - /GUI/			- GUI files and figures.
	- REST.m			- Contains many functions for necessary for the REST GUI
	- REST.fig		- Describes the REST GUI layout.
- /data/			- Demo datasets, channel locations, head models, and config files.
 - /dependencies/		- Necessary functions from MoBILAB and other sources. 

### References
L. Pion-Tonachini, S.-H. Hsu, S. Makeig, T.-P. Jung, G. Cauwenberghs, "Real-Time EEG Source-Mapping Toolbox (REST): Online ICA and Source Localization," in *IEEE EMBS*, 2015.

2015/08/26 Sheng-Hsiou Hsu / Luca Pion-Tonachini @ SCCN, UCSD
