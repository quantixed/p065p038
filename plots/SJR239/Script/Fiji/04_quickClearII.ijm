/*
 * Macro to use the glasbey view of all clusters segmented in 3D (and the log file of all images used)
 * to clear regions outside that are not of interest.
 */

glsby = getTitle();
no = getSliceNumber();
getDimensions(ww, hh, channels, slices, frames);
// add roi to roi manager
selectWindow(glsby);
roiManager("Add");
// hard-coded - sorry
Table.open("/Users/mlsmaf/Desktop/Log.txt");
// the log.txt file must have fpath at the beginning, 0-based
imgPath = Table.getString("fpath", no - 1);
// the text file comes from whatDoWeHave.ijm so the peaks_clust files are listed
// we want peaks_mod (the ones that come after nucleus deletion
imgPath = replace(imgPath, "_peaks_clust_", "_peaks_mod_");
//print(imgPath);
close("Log.txt");
// now open the corresponding image
open(imgPath);
seed = getTitle();
getDimensions(width, height, channels, slices, frames);
// calculate offset;
woff = (ww - width) / 2;
hoff = (hh - height) / 2;
// select roi in the image
roiManager("Select",0);
roiManager("translate", -woff, -hoff);
roiManager("Select",0);
// do the clearance
run("Clear Outside", "stack");
run("Select None");
save(replace(imgPath,"_peaks_mod_","_peaks_alt_"));
close();

//// get rid of roi
//roiManager("Delete");
//roiManager("update");
// empty the ROI manager
roiManager("reset");