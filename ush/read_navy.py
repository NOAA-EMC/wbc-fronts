import numpy
import glob
from datetime import datetime
import pandas as pd
#import pdb

# DCOM location
DCOMdir='/scratch2/NCEPDEV/marine/Todd.Spindler/noscrub/DCOM'

#dataDir='/marine/noscrub/Todd.Spindler/Navy'
hemisphere={'N':1,'S':-1,'E':1,'W':-1}

#----------------------------------------------------------------------------
def parse_contents(front,header,contents):
    front[header]={'lat':numpy.array([]),'lon':numpy.array([])}

    for line in contents:
        lat,lon=line[:5],line[5:]
        front[header]['lat']=numpy.append(front[header]['lat'],float(lat[:-1])*hemisphere[lat[-1]])
        front[header]['lon']=numpy.append(front[header]['lon'],float(lon[:-1])*hemisphere[lon[-1]])
    return front

#----------------------------------------------------------------------------
def read_navo(year=f'{datetime.now():%Y}'):
    """
    Decoder for the NAVO Gulf Stream Frontal Analysis text messages in /dcom
    """
    #pdb.set_trace()
    northname='GULF STREAM NORTH WALL'
    southname='GULF STREAM SOUTH WALL'
    eomname='FRONTAL DATA BASED ON MAX'
    
    front={}

    #year=datetime.now().strftime('%Y')
    # read in header and body
    #allfiles=glob.glob('/dcom/us007003/'+year+'*/wtxtbul/gs*.sub')
    allfiles=glob.glob(f'{DCOMdir}/{year}*/wtxtbul/gs*.sub')
    allfiles.sort()
    
    # check the date.  NAVO seems to be a day ahead.
    latest=allfiles[-1]
    latest_date=datetime.strptime(latest.split('/')[7],'%Y%m%d')

    for nf,filename in enumerate(allfiles):
        with open(filename) as f:
            contents=f.read().replace('\n','/').replace(':',' ').split('/')
            
            # remove blank lines
            contents=[token.strip() for token in contents if len(token.strip())>0]
            
            # split message into northwall and southwall parts    
            if ' '.join(contents).find(eomname) >=0:                
                [eom]=[nline for nline,line in enumerate(contents) if line.find(eomname) >=0]
            if ' '.join(contents).find(northname) >=0:
                [north]=[nline for nline,line in enumerate(contents) if line.find(northname) >=0]
            else:
                north=-999
            if ' '.join(contents).find(southname) >=0:
                [south]=[nline for nline,line in enumerate(contents) if line.find(southname) >=0]
            else:
                south=-999
            # sometimes north is first, sometimes south is first, sometimes south or north is missing
            # sometimes there are two norths or two souths.  Need to fix that rare case.
            northwall=[]
            southwall=[]
            if north <= south:
                if north != -999 and southwall != -999:
                    northwall=contents[north+1:south]
                    southwall=contents[south+1:eom]
                elif south==-999:
                    northwall=contents[north+1:eom]
                elif north==-999:
                    southwall=contents[south+1:eom]
            else:
                southwall=contents[south+1:north]
                northwall=contents[north+1:eom]
            thedate=datetime.strptime(''.join(contents[north].split()[-3:]),'%d%b%y')
            thedate=min(latest_date,thedate) # in case there's an offset in the date    
            if northname not in front:
                front[northname]={thedate:{}}
            if southname not in front:
                front[southname]={thedate:{}}
            northwall=' '.join(northwall).strip().split()
            southwall=' '.join(southwall).strip().split()
            front[northname]=parse_contents(front[northname],thedate,northwall)
            front[southname]=parse_contents(front[southname],thedate,southwall)
    
    #dict_of_df = {k: pd.DataFrame(v) for k,v in front.items()}
    #df = pd.concat(dict_of_df, axis=1)
    
    return front

#----------------------------------------------------------------------------
def read_navoceano(year=f'{datetime.now():%Y}'):
    """
    Decoder for the NAVOCEANO Ocean Frontal Analysis text message
    """
    
    front={}

    #year=datetime.now().strftime('%Y')
    # read in header and body
    #allfiles=glob.glob('/dcom/us007003/'+year+'*/wtxtbul/*xx.mrf')
    allfiles=glob.glob(f'{DCOMdir}/{year}*/wtxtbul/*xx.mrf')
    allfiles.sort()
    for filename in allfiles:
        with open(filename) as f:
            contents=f.read().replace('\n','/').split('/')
            
        # remove blank lines
        contents=[token.strip() for token in contents if len(token.strip())>0]
        
        # split message into header and body parts 
        header,body=contents[:21],contents[21:]
        
        # basic header syntax check
        if header[0] != 'APPROVED FOR PUBLIC RELEASE' or \
            header[1] != 'UNCLAS' or \
            header[3:6] != ['MSGID', 'NAVOCEANO', 'OVLY2'] or \
            header[8] != 'OVLY' or \
            header[13] != 'OCEAN FEATURE ANALYSIS':
            print('Header syntax error in '+filename)
            return
        
        # get the date from the header
        thedate=datetime.strptime(header[-1].split()[-1],'%d%b%y')
        
        # begin parse cycle for body (TEXT + LINE + LINE)
        while len(body) > 0 and body[0] != 'ENDAT':
            if body[0]=='ARC':
                body=body[7:]
                continue
            if body[0]=='TEXT':
                name,body=body[6],body[7:]
                if name not in front:
                    front[name]={}
                if body[0]=='LINE':
                    front[name][thedate]={'lat':numpy.array([]),'lon':numpy.array([])}
                continue
            if body[0]=='LINE':
                front[name][thedate]['lat']=numpy.append(front[name][thedate]['lat'],numpy.nan)  # put a break between lines
                front[name][thedate]['lon']=numpy.append(front[name][thedate]['lon'],numpy.nan)  # put a break between lines
                npoints,body=int(body[1]),body[4:]
                for np in range(npoints):
                    lats,lons,body=body[0],body[1],body[2:]
                    lat=float(lats[:4])/100.*hemisphere[lats[4]]
                    lon=float(lons[:5])/100.*hemisphere[lons[5]]
                    if lon < -180:
                        lon=lon%360
                    front[name][thedate]['lat']=numpy.append(front[name][thedate]['lat'],lat)
                    front[name][thedate]['lon']=numpy.append(front[name][thedate]['lon'],lon)
            else:
                break
            
    return front

#----------------------------------------------------------------------------
if __name__ == '__main__':
	
    n1=read_navo()
    n2=read_navoceano()
            
