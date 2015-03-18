function nut_checkver
% check for updates to NUTMEG

try
    currentver = urlread('http://nutmeg.berkeley.edu/currentversion.txt');
    currentver = currentver(1:(end-1)); % chop last character (should be carriage-return char);
catch
    errordlg('Sorry, update check failed. Something''s jacked with your Internet connection.');
    return;
end

if(strcmp(currentver,nut_ver))
    msgbox('NUTMEG is up-to-date.');
else
    msgbox(['Get with the times! NUTMEG version ' currentver ' is available from http://nutmeg.berkeley.edu']);
end

