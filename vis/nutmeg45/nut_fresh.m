% A script that clears all variables and functions, 
% closes all figure windows and cleans the console.

clear all;
try
    close all hidden;
end
clc;
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
nut_close;
