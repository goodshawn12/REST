REST (Real-time EEG Source-mapping Toolbox): Onlince ICA and Source Localization
--------------------------------------------------------------------------------

REST provides an easy-to-use and interactive platform to perform source-space analysis on live-streamed electroencephalographic (EEG) data in near real-time with visualizations of estimated independent source activations, corresponding source scalp topographies, source locations, and as well as other measures.


Files and Sub-directories:

 test_script.m		- Example script that shows how to run REST
 make_headModel.m	- Make your own head model for source localization
 bcilab_path.txt	- Define path to bcilab
 /functions		- All REST functions including online ICA and source localization
 /GUI			- GUI files and figures.
 /data			- Demo datasets, channel locations, head models, and config files.
 /dependencies		- Necessary functions from MoBILAB for source localizations. 


To use REST:

1. Make sure the BCILAB is called before running the main script. The latest versions of BCILAB can be found at github.com/sccn/BCILAB. When downloaded, go to bcilab directory and type 'bcilab' to initiate it.

2. Create LSL stream. LSL can be downloaded at ...


References:

L. Pion-Tonachini, S.-H. Hsu, S. Makeig, T.-P. Jung, G. Cauwenberghs, "Real-Time EEG Source-Mapping Toolbox (REST): Online ICA and Source Localization," in IEEE EMBS, 2015.




2014/11/06 Sheng-Hsiou Hsu / Luca Pion-Tonachini @ SCCN, UCSD

----------------

-----------
Instruction
-----------

