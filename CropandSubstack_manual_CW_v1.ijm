// This script helps to make cropping and substack generation a bit easier
// You still have to draw a box to define ROI to crop
//  And to specify start and end of substack

// this is the slice within which you want to generate the substack for 
// Needs to be manually changed
start =1;
end = 20;

// target directory also needs to be manually specified
target_dir = ""

run("Crop");
run("Make Substack...", "channels=1-3 slices=" + start + "-" + end);

Imagename = getTitle();

saveAs("Tiff", target_dir + Imagename);

close();
close();
