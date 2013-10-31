import numpy as np

def write1dmfarray(f, arr, nvalperline, textid):
    '''
    Write the modflow header and array to an open filed object.
    '''
    print '   writing array: ', textid
    f.write('internal 1 (free) -1 ' + textid + '\n')
    write1darray(f, arr, nvalperline)
    return

def write1darray(f, arr, nvalperline):
    '''
    Write the 1d array (arr) to the opened file.
    '''
    if arr.dtype == np.float: 
        sfmt = '{0:G} '
    elif arr.dtype == np.int32 or arr.dtype == 'int64':
        sfmt = '{0:n} '
    else:
        print 'type: ',arr.dtype
        raise Exception('unknown type in writearray')
    s = ''
    ivalonline = 0
    for i in xrange(len(arr)):
        val = arr[i]
        s += sfmt.format(val)
        ivalonline += 1
        if ivalonline == nvalperline or i == len(arr) - 1:
            s += '\n'
            f.write(s)
            s = ''
            ivalonline = 0
    return
    
def writejaarray(f, ia, jaarray, textid, nzero=0):
    '''
    Write the specified array to an open file object by
    putting one row on each line.  If the jaarray is
    the actual ja array, then set nzero to 1, and a value
    of one will be added to all values in the array.
    '''
    print '   writing array: ', textid
    f.write('internal 1 (free) -1 ' + textid + '\n')    
    if jaarray.dtype == np.float:
        sfmt = '{0:G} '
    elif jaarray.dtype == np.int:
        sfmt = '{0:n} '
    for n in xrange(len(ia) - 1):
        s = ''
        for m in xrange(ia[n], ia[n+1]):
            s += sfmt.format(jaarray[m] + nzero)
        f.write(s + '\n')
    return

def writebas(filename, qp, ibound, strt, hnoflo=-999.99, options=['free',]):
    '''
    Write the mfusg basic file
    '''
    print 'writing the bas file...'
    #open file
    f = open(filename, 'w')
    nlay = qp.modflowgrid.nlay
    
    #write options
    for opt in options:
        f.write(opt+' ')
    f.write('\n')
    
    #IBOUND(NDSLAY)
    istart = 0
    for k in xrange(nlay):
        istop = istart + qp.nodelay[k]
        write1dmfarray(f, ibound[istart: istop], 10, 'IBOUND(NDSLAY) K = ' + str(k + 1))
        istart = istop

    #write hnoflo
    f.write(str(hnoflo)+'\n')

    #STRT(NDSLAY)
    istart = 0
    for k in xrange(nlay):
        istop = istart + qp.nodelay[k]
        write1dmfarray(f, strt[istart: istop], 10, 'STRT(NDSLAY) K = ' + str(k + 1))
        istart = istop

    
def writedisu(qp, csrdata, filename, nper=1, stressperiods=[[1, 1, 1.0, 'SS']],
              top=None, bot=None):
    '''
    Write the unstructured discretization file
    '''
    print 'writing the disu file...'
    #open file
    f = open(filename, 'w')

    #NODES NLAY NJAG IVSD NPER ITMUNI LENUNI IDSYMRD
    nodes = csrdata.nodes
    nlay = qp.modflowgrid.nlay
    njag = csrdata.nja
    ivsd = 0
    itmuni = 1
    lenuni = 0
    idsymrd = 0
    
    s = '{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}{5:>10}{6:>10}{7:>10}\n'.format(
        nodes, nlay, njag, ivsd, nper, itmuni, lenuni, idsymrd)
    f.write(s)

    #LAYCBD(NLAY)
    s = ''
    for k in range(nlay):
        s += '0 '
    s += 'LAYCBD(NLAY)\n'
    f.write(s)
    
    #NODELAY(NLAY)
    write1dmfarray(f, qp.nodelay, 10, 'NODELAY(NLAY)')
    
    #Top(NDSLAY)
    if top is None:
        top = qp.top
    istart = 0
    for k in xrange(nlay):
        istop = istart + qp.nodelay[k]
        write1dmfarray(f, top[istart: istop], 10, 'TOP(NDSLAY) K = ' + str(k + 1))
        istart = istop
    
    #Bot(NDSLAY)
    if bot is None:
        bot = qp.bot
    istart = 0
    for k in xrange(nlay):
        istop = istart + qp.nodelay[k]
        write1dmfarray(f, bot[istart: istop], 10, 'BOT(NDSLAY) K = ' + str(k + 1))
        istart = istop
    
    #Area(NDSLAY)
    istart = 0
    for k in xrange(nlay):
        istop = istart + qp.nodelay[k]
        write1dmfarray(f, csrdata.areas[istart: istop], 10, 'AREA(NDSLAY) K = ' + str(k + 1))
        istart = istop
    
    #IAC(NODES)
    iac = np.empty( (nodes), dtype=int)
    for i in xrange(nodes):
        iac[i] = csrdata.ia[i + 1] - csrdata.ia[i]
    write1dmfarray(f, iac, 50, 'IAC(NODES)')
    
    #JA(NJAG)
    writejaarray(f, csrdata.ia, csrdata.ja, 'JA(NJAG)', nzero=1)
    
    #IVC(NJAG)
    #not writing because ivsd is not 1
    
    #CL12(NJAG)
    writejaarray(f, csrdata.ia, csrdata.cl1, 'CL12(NJAG)')
    
    #FAHL(NJAG)
    writejaarray(f, csrdata.ia, csrdata.fahl, 'FAHL(NJAG)')
    
    #PERLEN NSTP TSMULT SS/TR
    for l in stressperiods:
        p, n, t, s = l
        s = '{0} {1} {2} {3}\n'.format(p,n,t,s)
        f.write(s)
    f.close()
    return


def writelpf(filename, qp, ilpfcb=50, hdry=-888.88, nplpf=0, 
             ikcflag=0, options=['novfc', 'constantcv',], laytyp=None, 
             layavg=None, chani=None, layvka=None, laywet=None,
             wetfct=1.0, iwetit=1, ihdwet=0., anglex=None, hk=1.,
             hani=None, vka=1., ss=None, sy=None):
    '''
    Write the layer property flow input file
    '''
    print 'writing the lpf file...'
    #open file
    f = open(filename, 'w')
    nlay=qp.modflowgrid.nlay

    #ILPFCB HDRY NPLPF IKCFLAG [Options]
    line = str(ilpfcb) + ' '
    line += str(hdry) + ' '
    line += str(nplpf) + ' '
    line += str(ikcflag) + ' '
    for opt in options:
        line += opt + ' '
    line += '\n'
    f.write(line)
    
    #2. LAYTYP(NLAY)
    if laytyp is None: laytyp = np.zeros((nlay),dtype=np.int)
    write1darray(f, laytyp, 50)

    #3. LAYAVG(NLAY)
    if layavg is None: layavg = np.zeros((nlay),dtype=np.int)
    write1darray(f, layavg, 50)
    
    #4. CHANI(NLAY)
    if chani is None: chani = np.ones((nlay),dtype=np.int)
    write1darray(f, chani, 50)
    
    #5. LAYVKA(NLAY)
    if layvka is None: layvka = np.zeros((nlay),dtype=np.int)
    write1darray(f, layvka, 50)
    
    #6. LAYWET(NLAY)
    if laywet is None: laywet = np.zeros((nlay),dtype=np.int)
    write1darray(f, laywet, 50)
    
    #7. [WETFCT IWETIT IHDWET]
    #(Include item 7 only if LAYWET indicates at least one wettable layer.)
    if laywet.sum() > 0:
        line = str(wetfct) + ' ' + str(iwetit) + ' ' + str(ihdwet) + '\n'
        f.write(line)

    #8. [PARNAM PARTYP Parval NCLU] not supported
    
    #9. [Layer Mltarr Zonarr IZ] not supported

    #10. [ANGLEX(NJAG)] U1DREL. not supported

    istart = 0
    for k in xrange(nlay):
        istop = istart + qp.nodelay[k]
    
    #11. HK(NDSLAY) If any HK parameters are included, read only a print code.
        if np.isscalar(hk):
            f.write('constant ' + str(hk) + ' hk(NDSLAY) K = ' + str(k + 1) + '\n')
        else:
            write1dmfarray(f, hk[istart: istop], 10, 'hk(NDSLAY) K = ' + str(k + 1))

    #12. [HANI(NDSLAY)] Include item 12 only if CHANI is less than or equal to 0.
    #not supported
    
    #13. VKA(NDSLAY) If any VK or VANI parameters are included, read only a print code.
        if np.isscalar(vka):
            f.write('constant ' + str(vka) + ' vka(NDSLAY) K = ' + str(k + 1) + '\n')
        else:
            write1dmfarray(f, vka[istart: istop], 10, 'vka(NDSLAY) K = ' + str(k + 1))
    
    #14. [Ss(NDSLAY)] Include item 14 only if at least one stress period is transient. 
    #If there are any SS parameters, read only a print code
        if ss is not None:
            if np.isscalar(ss):
                f.write('constant ' + str(ss) + ' ss(NDSLAY) K = ' + str(k + 1) + '\n')
            else:
                write1dmfarray(f, ss[istart: istop], 10, 'ss(NDSLAY) K = ' + str(k + 1))

    #15. [Sy(NDSLAY)] Include item 15 only if at least one stress period is transient and 
    #LAYTYP is not 0. If any SY parameters are included, read only a print code.
    #enter sy for entire grid, even if not needed for some layers
        if sy is not None:
            if np.isscalar(sy):
                f.write('constant ' + str(sy) + ' sy(NDSLAY) K = ' + str(k + 1) + '\n')
            else:
                write1dmfarray(f, sy[istart: istop], 10, 'sy(NDSLAY) K = ' + str(k + 1))
    
    #16. [VKCB(NDSLAY)] Include item 16 only if IKCFLAG=0 and LAYCBD (in the Discretization 
    #File) is not 0. If any VKCB parameters are included, read only a print code.
    #not supported

    #17. [WETDRY(NDSLAY)] Include item 17 only if LAYWET is not 0 and LAYTYP is not 0 or 4.
    #not supported
    
    #If IKCFLAG is 1 or -1, indicating input of hydraulic conductivity (or transmissivity if 
    #confined) or inter-block conductance along connections then read item 18 for all connections 
    #over all layers. Otherwise, Item 18 is not read.
    #18. [Ksat(NJA)] U1DREL
    #not supported

        istart = istop
        
    return

def writewel(filename, iwelcb=50, wellist=[]):
    '''
    Write the well input file
    '''
    print 'writing the wel file...'
    #open file
    f = open(filename, 'w')

    #find max # of wells
    nper = len(wellist)
    mxactd = 0
    for wels in wellist:
        nwel = len(wels)
        if nwel > mxactd: mxactd = nwel
    
    #dataset 2
    line = '{0} {1}\n'.format(mxactd,iwelcb)
    f.write(line)
    
    for kper in xrange(nper):
        wels = wellist[kper]
        nwel = len(wels)
        line = '{0} {1}\n'.format(nwel, 0)
        f.write(line)
        for n, q in wels:
            line = '{0} {1}\n'.format(n, q)
            f.write(line)
    return
    
def writedrn(filename, idrncb=50, drnlist=[]):
    '''
    Write the drain input file
    '''
    print 'writing the drn file...'
    #open file
    f = open(filename, 'w')

    #find max # of drains
    nper = len(drnlist)
    mxactd = 0
    for drns in drnlist:
        ndrn = len(drns)
        if ndrn > mxactd: mxactd = ndrn
    
    #dataset 2
    line = '{0} {1}\n'.format(mxactd,idrncb)
    f.write(line)
    
    for kper in xrange(nper):
        drns = drnlist[kper]
        ndrn = len(drns)
        line = '{0} {1}\n'.format(ndrn, 0)
        f.write(line)
        for n, e, c in drns:
            line = '{0} {1} {2}\n'.format(n, e, c)
            f.write(line)
    f.close()
    return

def writeghb(filename, ighbcb=50, ghblist=[]):
    '''
    Write the ghb input file
    '''
    print 'writing the ghb file...'
    #open file
    f = open(filename, 'w')

    #find max # of ghbs
    nper = len(ghblist)
    mxactd = 0
    for ghbs in ghblist:
        nghb = len(ghbs)
        if nghb > mxactd: mxactd = nghb
    
    #dataset 2
    line = '{0} {1}\n'.format(mxactd,ighbcb)
    f.write(line)
    
    for kper in xrange(nper):
        ghbs = ghblist[kper]
        nghb = len(ghbs)
        line = '{0} {1}\n'.format(nghb, 0)
        f.write(line)
        for n, s, c in ghbs:
            line = '{0} {1} {2}\n'.format(n, s, c)
            f.write(line)
    f.close()
    return

def writerch(filename, irchcb=50, nrchop=1, rchlist=[]):
    '''
    Write the rch input file
    rchlist = [ [rech_kper1], [None], [rech_kper3], ...] where None will repeat
        recharge from last stress period, or
    rchlist = [ [rech_kper1, irch_kper1], [None, None], [rech_kper3, irch_kper3], ...]
        this is used for nrchop = 2.
    '''
    print 'writing the rch file...'
    #open file
    f = open(filename, 'w')

    #find max # of recharge records
    nper = len(rchlist)
    
    #dataset 2
    line = '{0} {1}\n'.format(nrchop, irchcb)
    f.write(line)
    
    for kper in xrange(nper):
        
        #get rech and irch out of rchlist and set inrch and inirch
        rechrecord = rchlist[kper]
        rech = rechrecord[0] #recharge is the first record in rchlist
        if rech is None:
            inrch = 0
        else:
            inrch = 1

        irch = None
        if nrchop == 2:
            irch = rechrecord[1] #recharge is the first record in rchlist
        inirch = 0
        if nrchop == 2:
            if irch is None:
                inirch = 0
            else:
                inirch = 1
        line = '{0} {1}\n'.format(inrch, inirch)
        f.write(line)

        if inrch != 0:
            if np.isscalar(rech):
                f.write('constant ' + str(rech) + ' recharge kper ' + 
                        str(kper + 1) + '\n')
            else:
                write1dmfarray(f, rech[:], 10, 'rech kper' + str(kper + 1))

        if nrchop == 2:
            if inirch != 0:
                if np.isscalar(irch):
                    f.write('constant ' + str(irch) + ' irch kper ' + 
                        str(kper + 1) + '\n')
                else:
                    write1dmfarray(f, irch[:], 10, 'irch kper' + str(kper + 1))
    f.close()
    return

def writegnc(filename, gnclist=[], i2kn=0, isymgncn=0):
    '''
    Write the gnc input file
    '''
    print 'writing the gnc file...'
    #open file
    f = open(filename, 'w')

    #dataset 1
    npgncn = 0
    mxgnn = 0
    ifalphan = 0
    ngncnpn = len(gnclist)
    #find the maximum number of adjacent alpha j
    mxadjn = 0
    for gnc in gnclist:
        n, m, alphajlist, fraclist = gnc
        if len(alphajlist) > mxadjn:
            mxadjn = len(alphajlist)        
    line = '{0} {1} {2} {3} {4} {5} {6}\n'.format(
            npgncn, mxgnn, ngncnpn, mxadjn, i2kn, isymgncn, ifalphan)
    f.write(line)
    
    #dataset 3
    for gnc in gnclist:
        n, m, alphajlist, fraclist = gnc
        ntoadd = mxadjn - len(alphajlist)
        if ntoadd > 0: #need to add more alphajs to bring up to mxadjn
            for i in xrange(ntoadd):
                alphajlist.append(m)
                fraclist.append(0.0)
            
        line = '{0} {1} '.format(n + 1, m + 1)
        for aj in alphajlist:
            line += '{0} '.format(aj + 1)
        for fj in fraclist:
            line += '{0} '.format(fj)
        line += '\n'
        f.write(line)
    f.close()
    return
