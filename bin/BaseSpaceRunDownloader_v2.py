from urllib.request import Request, urlopen, URLError
import json
import math
import sys
import os
import socket
import optparse


def arg_parser():
    parser = optparse.OptionParser()
    parser.add_option('-r', dest='runid', help='Run ID: required')
    parser.add_option('-a', dest='accesstoken', help='Access Token: required')
    (options, args) = parser.parse_args()
    try:
        if options.runid is None:
            raise Exception
        if options.accesstoken is None:
            raise Exception
    except Exception:
        print("Usage: BaseSpaceRunDownloader_vN.py -r <RunID>"
              " -a <AccessToken>")
        sys.exit()
    return options


def restrequest(rawrequest):
    request = Request(rawrequest)
    try:
        response = urlopen(request)
        json_string = response.read()
        json_obj = json.loads(json_string)
    except URLError as e:
        print('Got an error code:', e)
        sys.exit()
    return json_obj


def downloadrestrequest(rawrequest, path):
    dirname = RunID + os.sep + os.path.dirname(path)
    if dirname != '':
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
    request = (rawrequest)
    outfile = open(RunID + os.sep + path, 'wb')
    try:
        response = urlopen(request, timeout=1)
        outfile.write(response.read())
        outfile.close()
    except URLError as e:
        print('Got an error code:', e)
        outfile.close()
        downloadrestrequest(rawrequest, path)
    except socket.error:
        print('Got a socket error: retrying')
        outfile.close()
        downloadrestrequest(rawrequest, path)


options = arg_parser()
RunID = options.runid
AccessToken = options.accesstoken
request = ('http://api.basespace.illumina.com/v1pre3/runs/{}/files?'
           'access_token={}').format(RunID, AccessToken)
json_obj = restrequest(request)
totalCount = json_obj['Response']['TotalCount']
noffsets = int(math.ceil(float(totalCount)/1000.0))
hreflist = []
pathlist = []
filenamelist = []

for index in range(noffsets):
    offset = 1000*index
    request = (
        'http://api.basespace.illumina.com/v1pre3/runs/{}/files'
        '?access_token={}&limit=1000&Offset={}'
        ).format(RunID, AccessToken, offset)
    json_obj = restrequest(request)
    nfiles = len(json_obj['Response']['Items'])
    for fileindex in range(nfiles):
        href = json_obj['Response']['Items'][fileindex]['Href']
        hreflist.append(href)
        path = json_obj['Response']['Items'][fileindex]['Path']
        pathlist.append(path)

for index in range(len(hreflist)):
    request = (
        'http://api.basespace.illumina.com/{}/content?access_token={}'
    ).format(hreflist[index], AccessToken)
    print('downloading {}'.format(pathlist[index]))
    downloadrestrequest(request, pathlist[index])
