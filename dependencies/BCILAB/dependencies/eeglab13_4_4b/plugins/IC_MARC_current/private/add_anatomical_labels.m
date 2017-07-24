function eeg = add_anatomical_labels(eeg)
% add_anatomical_labels returns the EEGLab data structure that was given as input
% labels of anatomical areas added to the dipole fit of each IC.
%
% Input:
% EEGLab data structure with dipole fits for all ICs in the field dipfit.
%
% Output:
% eeg: The same EEGLab data structure that was given as input, but with
% labels of anatomical areas from the Talairach atlas added to the dipole
% fits of ICs.

% This function was built by Laura Froelich (laura.frolich@gmail.com)
% from part of set_fit_dipoles.m in BCILab.
% Author of set_fit_dipoles.m: Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

db = org.talairach.Database;
mfile_location=mfilename('fullpath');
talairach_location = [mfile_location(1:find(mfile_location == '/', 1, 'last')) 'talairach.nii'];
db.load(talairach_location);
for k=1:length(eeg.dipfit.model)
    try
        if(isempty(eeg.dipfit.model(1,k).posxyz))
            eeg.dipfit.label(1,k) = cell(UNKNOWN);
        else
            label = db.findLabel( eeg.dipfit.model(1,k).posxyz(1), ... % x
                eeg.dipfit.model(1,k).posxyz(2), ... % y
                eeg.dipfit.model(1,k).posxyz(3)  ... % z
                );
            if isempty(label)
                eeg.dipfit.label(1,k) = cell(UNKNOWN);
            else
                eeg.dipfit.label(1,k)  = cell(label.toString());
            end
        end
    catch
        eeg.dipfit.model(k).labels = {};
    end
end
