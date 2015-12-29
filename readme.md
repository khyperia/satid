# Satellite Identification

This project aims to complete the link in a chain of tools to identify satellites that appear in hobby astrophotographs.

![Example image with a satellite passing through it](http://i.imgur.com/pXHZ0kE.jpg)

## Building

The library [libpredict](https://github.com/la1k/libpredict) is the only dependency, other than the C standard library and math library. The source *should* be portable to all platforms, as I tried not to use any nonstandard functions (only a c99 compiler should be needed). If this is not the case, notify me and I can make it more portable.

    # Works on my system, yours may not. Basic CMakeLists.txt is also provided.
    gcc main.c -lpredict -lm -o satid

Quirks: Note that because of my portability efforts, the -time flag is very quirky. I believe it parses the time according to the current timezone, with the final number being an hour offset from the local timezone (yyyy-mm-dd hh:mm:ss tz). When printing, though, it uses UTC, so check to make sure that one makes sense.

## Usage

This is one program in a host of other programs that I use to identify satellites. Here's an example pipeline for an image I took a while back:

**Step 1** (using a built astrometry.net package): Plate-solve the image

    ./solve-field ~/sats/telescope.5-30.22-44-4.tiff

**Step 2** (using some image editor/etc): find the pixel coordinates of each end of the satellite path

    # 1392,-69
    # 0,303

**Step 3**: Convert pixels to RA/Dec

    ./wcs-xy2rd -w ~/sats/telescope.5-30.22-44-4.wcs -x 1392 -y -69.2503482919899
    # 186.7049466193,13.2306173047
    ./wcs-xy2rd -w ~/sats/telescope.5-30.22-44-4.wcs -x 0 -y 303.73492799999997
    # 186.5844359049,12.4408798662

**Step 4**: Download 3LE elements

I signed up for [space-track](https://www.space-track.org) and downloaded the 3le (note: not tle, the program expects a format where "the 0th line" is the name of the satellite) for a range of 1 day surrounding the target date. Note that the resulting data is over a megabyte in size.

**Step 5**: Obtain your observing coordinates/time

Google Maps is good for this, just copy the URL's lat/lon. Note it has to be pretty accurate - I think down to the arcminute, but I've never tested.

**Step 6**: Run this program on the RA/Dec pairs

    ./satid -lat 41.897717 -lon -85.8580815 -alt 300 -time "2014-05-30 22:44:04 -1" -tle ~/sats/telescope.5-30.22-44-4.tle.txt <<< "186.3398758412,12.4136527840 186.0392236285,13.1558875331"

This outputs all satellites within 1 degree (configuable by -angle flag) of that RA/Dec within -10 to 5 minutes of the specified time (the time is supposed to be end of exposure)

     -- RA/Dec -- 186.339876,12.413653
    Sat: 0 FENGYUN 1C DEB (30714)
    Angle: 0.394997
    Time: 2014-05-31 02:40:46:0.631402
    Sat: 0 DELTA 1 DEB (10706)
    Angle: 0.099283
    Time: 2014-05-31 02:40:46:0.175406
    Sat: 0 SL-24 R/B (28367)
    Angle: 0.135811
    Time: 2014-05-31 02:42:39:0.234039
     -- RA/Dec -- 186.039224,13.155888
    Sat: 0 FENGYUN 1C DEB (30714)
    Angle: 0.430584
    Time: 2014-05-31 02:40:48:0.585995
    Sat: 0 DELTA 1 DEB (10706)
    Angle: 0.146489
    Time: 2014-05-31 02:40:50:0.171490
    Sat: 0 SL-24 R/B (28367)
    Angle: 0.131117
    Time: 2014-05-31 02:42:42:0.042109

Note that there is some imprecision in the calculation for some reason. I believe this might be due to not accounting for atmospheric refraction. However, the matching satellite will be a consistent angle away - in this case, it's "SL-24 R/B (28367)". (the 0 in the text output is leftover from space-track's format - they include a 0 in it as 'the 0th line', other places don't. Just ignore it.) Also, order is consistent, in order of the tle file.

The result can be verified by opening the TLE file in [Stellarium](http://www.stellarium.org/) with the satellite plugin, and visually checking the angle, etc. - also useful is [calsky](http://www.calsky.com), where you can look up the brightness of a satellite and other information by the satellite number (28367 in this case). Finally, a sanity-check on the time is useful, to make sure it was actually within the exposure range.

Have fun doing astrophotography!
