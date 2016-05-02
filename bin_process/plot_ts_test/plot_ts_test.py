#!/usr/bin/env python2
"""
Tests plot_ts.sh by comparing known outputs (in TestCasePng) to ones generated.
Additional source file is passed to modify variables.
If there are visual differences:
    original images, inverted in area of change and diff image are shown.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 21 April 2016

USAGE: execute from current directory only './plot_ts_test.py' or 'python plot_ts_test.py'
    make sure python is version 2.x/2.6+

ISSUES: make executable from any location, could only show diff cropped and stretched.
"""

from shutil import copyfile
from subprocess import call
from os.path import basename, exists
from glob import glob
from PIL import Image

test_configs = [ \
'ts_start=0', \
'ts_inc=53', \
'ts_total=4']
params = 'test.par'
default_params = 'test_default.par'

copyfile(default_params, params)
fp = open(params, 'a')
fp.write('\n'.join(test_configs) + '\n')
fp.close()

def invert_diff_area(image, diff_box):
    if diff_box[0] > diff_box[1] or diff_box[2] > diff_box[3]:
        return

    for x in xrange(diff_box[0], diff_box[1] + 1):
        for y in xrange(diff_box[2], diff_box[3] + 1):
            pixel = image.getpixel((x, y))
            image.putpixel((x, y), (255 - pixel[0], 255 - pixel[1], 255 - pixel[2]))

def locate_diff_box(pl1, pl2, width, height):
    diff_box = [width, 0, height, 0]

    # find top
    scanner = 0
    while pl1[scanner] == pl2[scanner]:
        scanner += 1
    diff_box[2] = scanner / width

    # find bottom
    scanner = width * height - 1
    while pl1[scanner] == pl2[scanner]:
        scanner -= 1
    diff_box[3] = scanner / width

    # find left
    scanner = 0
    while pl1[scanner] == pl2[scanner]:
        scanner += width
        if scanner / width >= height:
            scanner -= width * height
            scanner += 1
    diff_box[0] = scanner % width

    # find right
    scanner = width * height - 1
    while pl1[scanner] == pl2[scanner]:
        scanner -= width
        if scanner < 0:
            scanner += width * height
            scanner -= 1
    diff_box[1] = scanner % width

    return diff_box

def highlight_change(img1, img2):
    width, height = img1.size

    img_diff = Image.new('RGB', (width, height))
    img1_highlight = Image.new('RGB', (width, height))
    img2_highlight = Image.new('RGB', (width, height))
    img_diff.paste(img1)
    img1_highlight.paste(img1)
    img2_highlight.paste(img2)

    pl1 = list(img1.getdata())
    pl2 = list(img2.getdata())
    diff_box = locate_diff_box(pl1, pl2, width, height)

    for p in xrange(diff_box[2] * width + diff_box[0], diff_box[3] * width + diff_box[1] + 1):
        if pl1[p] != pl2[p]:
            # arbitrary factors used to balance contrast/visibility
            diff = 255 - int(((abs(pl1[p][0] - pl2[p][0]) + \
                    abs(pl1[p][1] - pl2[p][1]) + \
                    abs(pl1[p][2] - pl2[p][2])) / 9) ** 2 * 0.6)
            if diff < 0:
                diff = 0
            diff += 102
            if pl1[p] > pl2[p]:
                pixel = (diff, 0, diff)
            else:
                pixel = (0, diff, diff)
            img_diff.putpixel((p % width, p / width), pixel)

    invert_diff_area(img1_highlight, diff_box)
    invert_diff_area(img2_highlight, diff_box)
    return img1_highlight, img2_highlight, img_diff

def combine_images(img1, img2, img1_h, img2_h, diff_h):
    width, height = img1.size
    width_thumb = 300
    height_thumb = int((float(width_thumb) / width) * height)

    img1.thumbnail((width_thumb, height_thumb))
    img2.thumbnail((width_thumb, height_thumb))
    img1_h.thumbnail((width_thumb, height_thumb))
    img2_h.thumbnail((width_thumb, height_thumb))
    diff_h.thumbnail((width_thumb * 2, height_thumb * 2))

    im_thumbs = Image.new('RGB', (width_thumb * 4, height_thumb * 2))
    im_thumbs.paste(img1, (0,0))
    im_thumbs.paste(img2, (width_thumb * 1, 0))
    im_thumbs.paste(img1_h, (0, height_thumb))
    im_thumbs.paste(img2_h, (width_thumb, height_thumb))
    im_thumbs.paste(diff_h, (width_thumb * 2, 0))
    return im_thumbs

def visualisation_manual():
    print('||=================================||')
    print('||          |         |  TESTCASE  ||')
    print('|| TESTCASE | RESULT  |     VS     ||')
    print('||          |         |   RESULT   ||')
    print('||---------------------  DIFF IMG  ||')
    print('|| TESTCASE | RESULT  |  ADD/DELL  ||')
    print('||   DIFF   |  DIFF   | A: MAGENTA ||')
    print('||   AREA   |  AREA   | D: CYAN    ||')
    print('||=================================||')




print('Running plot_ts on local test case...')
call(['bash', '../plot_ts.sh', '4', 'test.par'])

testcases = map(basename, glob('TestCasePng/*'))
manual_displayed = False
for case in testcases:
    if not exists('Png/' + case):
        print('WARNING: png file not created: ' + case)
        continue

    img1 = Image.open('TestCasePng/' + case)
    img2 = Image.open('Png/' + case)

    if img1 == img2:
        print('Success: Images are the same: ' + case)
        continue

    if img1.size != img2.size:
        print('WARNING: Image produced is in different resolution: ' + case)
        print('Will not produce diff for differently sized images.')
        continue

    print('WARNING: Image produced is different, producing diff for: ' + case)
    if not manual_displayed:
        visualisation_manual()
        manual_displayed = True
    print('Rendering...')

    img1_h, img2_h, diff_h = highlight_change(img1, img2)
    visualisation = combine_images(img1, img2, img1_h, img2_h, diff_h)
    visualisation.show()
    print('Rendering complete.')


print('Process Complete.')


