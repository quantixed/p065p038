/*
 * Routine to take composites generated from nnU-Net inference files (0,1,2)
 * make a large stack of the substacks and a nice stack to visualise
 */

dir = File.directory();
stackList = getList("image.titles");
filename = "stk_" + stackList[0]
if(stackList.length > 1) {
	run("Concatenate...", "all_open title=new");
}
save(dir + filename);
run("Z Project...", "projection=[Average Intensity]");
Stack.setChannel(1);
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(2);
run("Enhance Contrast", "saturated=0.35");
filename = "avg_" + stackList[0];
save(dir + filename);
run("RGB Color");
filename = replace(filename, ".tif", ".jpg");
saveAs("jpg", dir + filename);
close("*");