'''
This module will (hopefully) read MSA files and color regions on them according to a list
of coordinates entered for each sequence
'''
from PyRTF import *
from Bio import AlignIO

class hMSA(object):
    '''
    Holds a multiple sequence alignment object (Bio.Align.MultipleSeqAlignment)
    Receives as well regions in the different sequences (key is the id) to be 
    then colored in the pretty plotting of the MSA
    '''

    def __init__(self,msa):
        self.msa=msa
        self.regions=[]
        self.ids={}
        for aln in msa:
            self.ids[aln.id]=aln

        
    def addRegion(self,dict):
        '''
        Map the regions in the dictionary to their coordinates in the MSA.
        The value is a list of tupples, each with the start and end of the
        region.
        The values in the dictionary are the from and to coordinates python 
        style (the end is actually 1-based, not 0-based like the start)
        '''
        dict2={}
        for key,val in dict.items():
            dict2[key]=[]
            for pos1,pos2 in val:
                dict2[key].append((self.mapGaps(key,pos1),self.mapGaps(key,pos2-1)+1))
        self.regions.append(dict2)
        
    def locateID(self,key):
        if key in self.ids:
            return self.ids[key]

    def mapGaps(self,key,pos):
        aln = self.locateID(key)
        if aln:
            gpos=-1
            seqpos=0
            while gpos<pos:
                if aln[seqpos]!='-': gpos+=1
                seqpos+=1
            return seqpos-1

    def realPos(self,key,pos):
        '''
        Return the real position along a sequence (remove gaps)
        '''
        aln=self.locateID(key)
        if aln: return pos+1-aln.seq[:pos].count('-')

    def getPosColors(self,key,pos):
        '''
        return -1 if not colored, 0 for the first mark, 1 for the second and 2 for both 
        '''
        colorkey=0
        for i,colors in enumerate(self.regions):
            if key in colors:
                for pos1,pos2 in colors[key]:
                    if pos in range(pos1,pos2):
                        colorkey=(i+1)
        return colorkey

        
    def getPosNuc(self,key,pos):
        aln=self.locateID(key)
        if aln:
            return aln[pos]
        

def print2RTF(hMSA,file,txids=None,width=40):
    '''
    print the MSA to an RTF file with some regions in different colors, as defined in the hMSA
    object. The colors are set arbitrarily to redblue and yellow, with green in the intersection
    '''
    import txid2path
    if not txids:
        txids=hMSA.ids.keys()
    #preparing the document
    doc     = Document()
    ss      = doc.StyleSheet
    section = Section()
    doc.Sections.append( section )
    colors=[ss.Colours.White, ss.Colours.Yellow, ss.Colours.Blue, ss.Colours.Green, ss.Colours.Red, ss.Colours.Pink, ss.Colours.Violet]
    slen = len(hMSA.msa[0])
    spos = 0
    font=ss.Fonts.CourierNew
    while spos<=slen:
        for txid in txids:
            #print the header
            p=Paragraph()
            p.append(TEXT(txid2path.getInitials(txid),font=font),TAB,TEXT(str(hMSA.realPos(txid,spos)),font=font),TAB)
            for pos in range(spos,min(spos+width,slen)):
                p.append(Text(hMSA.getPosNuc(txid,pos),TextPS(font=font),ShadingPS(background=colors[hMSA.getPosColors(txid,pos)])))
            section.append(p)
        spos+=width
        #add an empty line
        p=Paragraph()
        section.append(p)
    DR = Renderer()
    DR.Write( doc, file  )


    

