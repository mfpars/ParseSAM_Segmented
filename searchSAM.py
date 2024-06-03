#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:36:55 2024

@author: mfpars
"""
import os
import pandas as pd
from pandas import json_normalize
import re


import matplotlib.pyplot as plt
import seaborn as sns 
from matplotlib.lines import Line2D

#%% Define functions

def searchSAM(base_dir,samfilename,reference_name, term='segmented'):

    sampath = os.path.join(base_dir,  samfilename)
    out_file = os.path.join(base_dir, samfilename + term +'_found_readnames.csv')

   
    # initiate output dataframe
    # found = pd.DataFrame(columns=['QNAME', # Query name of read or read pair
    #                               'FLAG', # bitwise FLAG, indicates pairing, strand, etc
    #                               'RNAME', # Reference sequence name
    #                               'POS', # 1-based leftmost position of clipped alignment
    #                               'MAPQ', # Phred-scaled Mapping quality
    #                               'CIGAR', # extended CIGAR string
    #                               'RNEXT', # Mate/next read reference name ('=' if same as RNAME)
    #                               'PNEXT', # 1-based leftmost position of mate/next read alignment
    #                               'ISIZE', # inferred insert size
    #                               'SEQ', # Query sequence on the same strand as the reference
    #                               'QUAL', # Phred base quality of Query
    #                               'Tags', # list of extra tags
    #                               ])
    ind = 0 # initiate output index
    found_entries = [] # initiate output index
    with open(sampath, 'r') as sam:
        for line in sam:
            if line.startswith('@'):
                continue  # Skip header lines
                
            line = line.strip() # remove newline characters and spaces from the end of the line
            fields = line.split('\t') # Split the line by tab symbol into individual fields
            flag = int(fields[1])
            
            # Filter for primary alignments (flag not containing bit 0x100) 
            # reads mapped to the desired reference (e.g. viral genome)
            # and whose associated next read segment also maps to the reference, if anything
            if (not flag & 0x100 and
                fields[2] == reference_name and
                (fields[6] == reference_name or fields[6] == '*')):
            
                # The read name is the first field
                # read_name = fields[0]
                
                if term == 'segmented':
                    # Check for the presence of "TC", "FC", or "FI" tags in the
                    # last set of fields; this indicates a chimeric read, possibly a copyback
                    for field in fields[11:]:
                        if field.startswith("TC") or field.startswith("FC") or field.startswith("FI"):
                            # # save read name to file
                            # with open(out_file, 'a') as output_file:
                            #     output_file.write(read_name + '\n')
                            # save SAM entry to copybacks dict
                            found_entries.append(fields[0:11]+[fields[11:]]) # group all TAG fields into single Tags column
                            
                            ind+=1 # iterate index
                            break
                        if ind % 100 == 0: # just to update that it is running
                            print(ind, ' segmented entries')
                            
    
    # note: much faster to append to list in for loop then convert to DataFrame than to directly append to DataFrame in for loop.
    found = pd.DataFrame(found_entries,
                         columns=['QNAME', # Query name of read or read pair
                                'FLAG', # bitwise FLAG, indicates pairing, strand, etc
                                'RNAME', # Reference sequence name
                                'POS', # 1-based leftmost position of clipped alignment
                                'MAPQ', # Phred-scaled Mapping quality
                                'CIGAR', # extended CIGAR string
                                'RNEXT', # Mate/next read reference name ('=' if same as RNAME)
                                'PNEXT', # 1-based leftmost position of mate/next read alignment
                                'ISIZE', # inferred insert size
                                'SEQ', # Query sequence on the same strand as the reference
                                'QUAL', # Phred base quality of Query
                                'Tags', # list of extra tags
                                ])
    found = found.astype({'FLAG': int, 'POS': int, 'PNEXT': int})
    return found


def copyback_filter(segments_df):
# Filter for copybacks specifically
    # Determine Reverse Complementarity from the FLAG
    segments_df['RevComp'] = segments_df['FLAG'] & 0x10

    # Group by QNAME and calculate the RevComp conditions
    # any is True if any RevComp is True, and all is True if all RevComp are True for each QNAME.
    group_stats = segments_df.groupby('QNAME')['RevComp'].agg(['any', 'all'])
    
    # Determine which QNAMEs have both RevComp True and False (i.e., are copybacks)
    # (any is True means at least one RevComp is True, and ~all is True means not all are True, so there's at least one False RevComp)
    copyback_qnames = group_stats[(group_stats['any']) & (~group_stats['all'])].index
    
    # Create a boolean Series indicating whether each row's QNAME is in the copyback_qnames
    copybacks_loc = segments_df['QNAME'].isin(copyback_qnames)
    
    # readnames = (q for q in segments_df['QNAME'].unique()) # generator for unique read names
    # 
    # segments_df['RevComp'] = segments_df['FLAG'].apply(lambda f: f & 0x10)
    # 
    # copyback_q = {}
    # for q in readnames:
    #     # Only if the segments for read q contain both RevComp and not Revcomp (True and False values),
    #         # there is a possible copyback.
    #     if not segments_df.loc[segments_df['QNAME'] == q]['RevComp'].all( # if not all of the RevComp values for q are True
    #             ) and not (~segments_df.loc[segments_df['QNAME']==q]['RevComp']).all(): # and not all of the RevComp values for q are False
    #         copyback_q[q] = True
    #     else:
    #         copyback_q[q] = False
    
    # copybacks_loc = segments_df['QNAME'].apply(lambda x: copyback_q[x])
    return copybacks_loc


    # Flag notes:
        # bit 16 (0x10) indicates whether SEQ is reverse complemented
        # bit 32 (0x20) indicates whether SEQ of the next segment in the read is reverse complemented


def melt_segments(SAM_df):
    # SAM_df is a dataframe in SAM format
    
    # Helper function to extract tag values
    def tagn(TAGS,searchtag):
        for item in TAGS:
            if item.startswith(searchtag):
                return int(item.split(':')[2]) # tag format is FI:i:n, where n is the segment id

    # Add the segment id (FI) and total segment count (TC) from the tags as new columns in the df
    SAM_df['FI'] = SAM_df['Tags'].apply(lambda tags: tagn(tags, 'FI'))
    SAM_df['TC'] = SAM_df['Tags'].apply(lambda tags: tagn(tags, 'TC'))

    # Initialize lists/columns of outputs
    # Initialize new columns
    # SAM_df['Start'] = 0
    # SAM_df['End'] = 0
    # SAM_df['Length'] = 0
    # SAM_df['RevComp'] = False
    # melted_reads = []
    # MINDs = [] # Mutation, Insertion, (N)recombination, Deletions

    ## Pre-compile regex expressions    
    
    # Regex pattern to split at the border between letters and numbers,
    # i.e. between each alignment section (soft padding, matching, deletions, etc)
    cigarsplit_regex = re.compile(r'(?<=[A-Za-z])(?=\d)')
    # Regex pattern to split at the border between numbers and letters,
    # i.e. to separate the number of bp and the alignment type
    chunksplit_regex = re.compile(r'(?<=\d)(?=[A-Za-z])')
    
    
    def parse_cigar(cigar):
        cigar_chunks = re.split(cigarsplit_regex, cigar)
        return [re.split(chunksplit_regex, chunk) for chunk in cigar_chunks]

    def segparse(chunked_cigar):
        # identify the start and end points of the aligned section of sequence
        seglength = 0
        start_idx = 0
        seg_mind_pos = []
        
        for idx, item in enumerate(chunked_cigar):
            event_length = int(item[0])
            if item[1] in ['M', 'X']:
                seglength += event_length
                start_idx = idx + 1
                break
        
        for item in chunked_cigar[start_idx:]:
            event_length = int(item[0])
            if item[1] == 'H': #hard clip, indicates end of segment alignment
                break
            if item[1] in ['X', 'D', 'N', 'M']: # mutations, deletions/splices, and matches continue the count along the reference genome
                if item[1] != 'M': # not an exact match, record MIND event
                    seg_mind_pos.append({'MINDtype': item[1],
                                         'Start': seglength, # start position of event relative to segment alignment start
                                         'Length': event_length,
                                         'End': seglength + event_length})
                seglength += event_length
            elif item[1] in ['I', 'S']:
                # insertions don't add to the count in the reference genome position
                # but still want to save this information.
                # Virema records large insertions as soft pads ('S') in the CIGAR string

                seg_mind_pos.append({'MINDtype': item[1],
                                     'Start': seglength,
                                     'Length': event_length,
                                     'End': seglength}) # relative to reference genome, insertion starts where it ends

        return seglength, seg_mind_pos

    # Vector operations on Pandas columns to calculate Start/End/Length/RevComp
    SAM_df['Start'] = SAM_df['POS'].astype(int)
    
    chunks = SAM_df['CIGAR'].apply(parse_cigar)
    segparse_results = chunks.apply(segparse)
    aligned_length, relative_mind_pos = zip(*segparse_results)
    
    SAM_df['End'] = SAM_df['Start'] + aligned_length
    SAM_df['Length'] = aligned_length
    SAM_df['RevComp'] = (SAM_df['FLAG'] & 0x10) != 0
    
    relative_minds_df = pd.DataFrame({'minds': relative_mind_pos})
    relative_minds_df['QNAME'] = SAM_df['QNAME'].values
    relative_minds_df['Seg'] = SAM_df['FI'].values
    relative_minds_df['POS'] = SAM_df['POS'].values
    
    # Explode the 'minds' column to separate each dictionary into its own row
    relative_minds_exploded = relative_minds_df.explode('minds')
    relative_minds_normal = json_normalize(relative_minds_exploded['minds'])
    
    MINDs = pd.concat([relative_minds_exploded[['QNAME', 'Seg', 'POS']].reset_index(drop=True),
                   relative_minds_normal], axis=1)
    
    MINDs['Start'] = MINDs['Start'] + MINDs['POS'] # update start to genome location
    MINDs['End'] = MINDs['End'] + MINDs['POS'] # update end to genome location
    MINDs.drop(columns='POS', inplace=True)
    return SAM_df, MINDs.dropna()

#%%
def plotreads(base_dir, samfilename, reference_name, parsed_reads_df, MINDs=[]):
    ### Plot set up
            
        ## Define a y-value to space out the reads on the plot
    
        read_ind = {qname: idx for idx, qname in enumerate(parsed_reads_df['QNAME'].unique())}
        parsed_reads_df['Read_Index'] = [read_ind[row['QNAME']] for idx, row in parsed_reads_df.iterrows()]
        
        unique_reads = parsed_reads_df.drop_duplicates('QNAME')[['QNAME','TC']]
        unique_reads['read_y'] = unique_reads['TC'].cumsum()
        read_y = {row['QNAME']: row['read_y'] for ind, row in unique_reads.iterrows()}
        parsed_reads_df['y'] = [4*read_y[row['QNAME']] + row['FI'] for idx, row in parsed_reads_df.iterrows()]
        
        # parsed_reads_df['y'] = 50*parsed_reads_df['Read_Index'] + 1*parsed_reads_df['FI']
        y_seg = {(row['QNAME'], row['FI']): row['y'] for idx, row in parsed_reads_df[['QNAME','FI','y']].iterrows()}
        # parsed_reads_df['y'] = [y_read[row['QNAME']] + 6*int(row['FI']) for idx, row in parsed_reads_df.iterrows()]
    
        ## Define figure parameters

        width = (parsed_reads_df['End'].max()) * 0.01 
        height = len(parsed_reads_df) * 0.05 # taller plot if more segments to plot
        fig, ax = plt.subplots(figsize=(width, height), layout="constrained")
        colors = sns.color_palette(# palette = 'gist_earth',
                                   n_colors = 5)
        strandcolor = {True: colors[0], False: 'black'} # set up dict to define color based on RevComp being true or false
    
    ###  Plot segments and connectors
    
        for index, row in parsed_reads_df.iterrows():
            ax.hlines(y=row['y'], xmin=row['Start'], 
                       xmax=row['End'], 
                       color=strandcolor[row['RevComp']],
                       linewidth = 0.5,
                       )
        
        for qname in parsed_reads_df['QNAME'].unique():
            
            
            # connect = [] # empty list to populate with segment connector
            for i in parsed_reads_df.loc[parsed_reads_df['QNAME']==qname]['FI']: 
                # TC is total number of segments
                # FI (index of segment in read) is 1-indexed
                
                
                segment = parsed_reads_df.loc[(parsed_reads_df['QNAME']==qname) &
                                              (parsed_reads_df['FI']==i)]
                
                nextsegment = parsed_reads_df.loc[(parsed_reads_df['QNAME']==qname) &
                                              (parsed_reads_df['FI']==i+1)]
                
                if nextsegment.empty: # reached the end of the read
                
                    # write read index on the plot
                    # ax.text(0, segment['y'], read_ind[segment['QNAME'].values[0]],
                    #         horizontalalignment='right',
                    #         verticalalignment='center',
                    #         **{'fontfamily': 'sans-serif','fontsize': 8})
                
                    continue
                
                # find start point of connector: 
                    # the End of a forward strand segment, or
                    # the Start of a reverse-complemented segment
                revcomp = segment['RevComp'].iloc[0]
                if revcomp:
                    connectstart = segment['Start'].iloc[0]
                else:
                    connectstart = segment['End'].iloc[0]
                
                # find end point of connector: 
                    # the Start of a forward-stranded next segment, or
                    # the End of a reverse-complemented next segment
                revcomp_next = nextsegment['RevComp'].iloc[0]
                if revcomp_next:
                    connectend = nextsegment['End'].iloc[0]
                else:
                    connectend = nextsegment['Start'].iloc[0]
                
                startpoint = {'QNAME': qname,
                              'pos': connectstart,
                              'y': segment['y']
                              }
                endpoint = {'QNAME': qname,
                            'pos': connectend,
                            'y': nextsegment['y']
                            }
                connector = pd.DataFrame([startpoint,endpoint])
                 
                ax.plot(connector['pos'],connector['y'],
                        # linestyle = '--',
                        dashes = [2,1],
                        linewidth = 0.25,
                        color = 'grey',
                        marker = "none"
                     )              
                # connect.append(startpoint)
                # connect.append(endpoint)
            
                
            
    ### Plot mutation/insertion/deletion events
    
        if not MINDs.empty:  # if events were identified
            MINDs['y'] = [y_seg[(row['QNAME'],row['Seg'])] for idx, row in MINDs.iterrows()]
            
            # Plot deletions as colored horizontal bars over the aligned segment
            Deletions = MINDs.loc[(MINDs['MINDtype']=='D') | (MINDs['MINDtype'] == 'N')]
            for index, row in Deletions.iterrows():
                ax.hlines(y=row['y'], xmin=row['Start'], xmax=row['End'], 
                           color=colors[1],
                           linewidth=1,
                           linestyle = ':',
                           gapcolor='grey',
                           alpha = 0.6,
                           # ax = ax
                           )
            
            # Plot Insertions and Mutations as point markers over the aligned segment bars
            MIN = MINDs.loc[(MINDs['MINDtype'] == 'X') | (MINDs['MINDtype'] == 'S')
                               | (MINDs['MINDtype'] == 'I')]
            
            sns.scatterplot(data = MIN,
                            x = 'Start',
                            y = 'y',
                            hue = 'MINDtype',
                            style = 'MINDtype',
                            hue_order = ['X','S','I'],
                            markers = {'X': 'X', 'S': 'o', 'I': 's'},
                            ax = ax,
                            palette = colors[2:],
                            s = 5, # adjust size of the MIND markers on the plot
                            alpha = 0.3, # make transparent to visualize lines underneath
                            )
            
    ### Beautify the plot for visual clarity and save figure
        
        legend_elements = [Line2D([0], [0], color=colors[0], lw=2, label='Reverse complement mapped'),
                           Line2D([0], [0], color='black', lw=2, label='Forward strand mapped'),
                           Line2D([0], [0], color='grey', lw=1, label='Segment connection',
                                  linestyle = '--'),
                           Line2D([0], [0], color=colors[1], lw=4, linestyle = ':',
                                  gapcolor='grey', alpha = 0.6, label='Deletion'),
                           Line2D([0], [0], marker='x', color=colors[2], label='Mutation',
                                  lw = 0, markersize=5, alpha = 0.3),
                           Line2D([0], [0], marker='o', color=colors[3], label='Soft pad',
                                  lw = 0, markersize=5, alpha = 0.3),
                           Line2D([0], [0], marker='s', color=colors[4], label='Insertion',
                                  lw = 0, markersize=5, alpha = 0.3),
                           ]
    
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1, 1, 0.25),
                  frameon = False,
                  handles=legend_elements
                  )
        ax.set_xlabel('Genome Position')
        ax.set_ylabel('Read')
        
        # Set y-ticks and labels
        first_seg_y = [parsed_reads_df[parsed_reads_df['QNAME'] == qname].iloc[0]['y'] for qname in parsed_reads_df['QNAME'].unique()]
        ax.set_yticks(first_seg_y, labels = [read_ind[qname] for qname in parsed_reads_df['QNAME'].unique()])
        ax.tick_params(width = 0.5, color = 'grey')
        ax.set_ylim(0,parsed_reads_df['y'].max())
        # ax.set_yticklabels(read_ind[qname] for qname in parsed_reads_df['QNAME'].unique())
        # ax.set_yticks([]) # remove arbitrary y-axis ticks and tick-labels, which are just counting arbitrary spacing on the plot
        ax.spines['bottom'].set_position('zero')
        ax.spines['left'].set_position('zero')
        sns.despine(ax=ax,top = True, right = True, left = True, bottom = False)

        
        # Add verticle grid lines at x-tick positions
        ax.grid(visible = True, axis = 'x', color='grey', linestyle='-', linewidth=2,
                alpha = 0.2)
        
        # Adjust subplot parameters to make room for labels
        # fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        
        outfolder = os.path.join(base_dir, samfilename[:-4])
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        fig.savefig(os.path.join(outfolder,reference_name + '_Copyback_diagram.pdf'),bbox_inches = 'tight'
                    )
        
        
        ### Save csv file to refer plot read index to read name
        
        parsed_reads_df.to_csv(os.path.join(outfolder,reference_name+'Segmented_reads_parsed.csv'))

#%%    
def calc_dsRNA(melted_reads_df, ):
    # Filter melted_reads dataframe for reads containing doublestranded regions
    # (i.e. copybacks), and calculate how many double-stranded base pairs the read has.
    
    readnames = (q for q in melted_reads_df['QNAME'].unique()) # collect unique readnames
    
    for q in readnames: # for each read
        read_rows = melted_reads_df.loc[melted_reads_df['QNAME'] == q]
        
        # if the read contains both forward and reverse segments (RevComp = True for one segment,
            # RevComp = False for another segment), there is a chance of being double-stranded
        if (True in read_rows['RevComp']) and (False in read_rows['RevComp']):
            TC = read_rows.loc[0,'TC'] # total number of segments in the read
            
            for seg in range(1,TC): # segment is 1-indexed, so start at 1
               
                # if this segment is in the dataframe and if there is a next segment:
                    # (the segment seg may map to a different genome as a result of recombination,
                    # and therefore not show up in this dataframe)
                if (seg in read_rows['Seg']) and (seg+1 in read_rows['Seg']):
                    segment = read_rows.loc[read_rows['Seg'] == seg]
                    next_segment = read_rows.loc[read_rows['Seg'] == seg+1]
                    
                    ''' Several possible conditions lead to hairpin potential (i.e. double-stranded base pairs) within the read:
                            1. Segment 1 is reverse-complemented and segment 2 is forward-stranded, and
                                a. The start of segment 2 is between start/end of segment 1, or
                                b. The start of segment 1 is between start/end of segment 2
                                
                            2. Segment 1 is forward-stranded and segment 2 is reverse-complemented, and
                                a. The start of segment 2 is between start/end of segment 1, or
                                b. The start of segment 1 is between start/end of segment 2
                        
                        (note, the 'start' and 'end' of segments listed here refers to the start/end of the matching aligned region,
                          so does not include soft pads or hard-clipped regions)
            
                    '''
                    if segment['RevComp'] != next_segment['RevComp']: # must map to opposite strands to be complementary
                       
                        if (segment['Start'] <= next_segment['Start'] <= segment['End']): # condition a.
                            
                        
                        
                        
                            #ds_region_start = min(segment['Start'], next_segment['Start']) # leftmost coordinate; depends on RevComp of segments
                            # ds_region_end = min()
                            # ds_region = {'Start': ds_region_start, 'ds_length (bp)': next_segment['Start']
                        
                        elif (next_segment['Start'] <= segment['Start'] <= next_segment['End'] # condition b.
                              ):

                        

#%% Identify segmented reads

if __name__ == '__main__':
    base_dir = '/Users/mfpars/RNAseq_data/fClip-opt_virema/20240313_Virema_wHost/dox1in'
    samfiles = ['dox1in_recombinations.SAM']
    references = ['HIV_genome_nonif']
    term = 'segmented'
    
    for file in samfiles:
         samfilename = file
         
         for reference_name in references:
             segmented = searchSAM(base_dir, samfilename, reference_name, 'segmented')
             segmented['Tags'] = segmented['Tags'].apply(tuple) # covert list of tags to hashable tuples
             
         # Due to how the filtering worked in above, the copybacks dataframe has
             # some lines that are the final segment of a read which initially partially mapped
             # to something other than the desired reference.
             # To further filter for only reads with multiple segments mapping to the specified reference,
             # select only rows whose "QNAME" value is duplicated elsewhere in the dataframe.
             # Identify duplicates in column 'A'
             duplicates_mask = segmented['QNAME'].duplicated(keep=False)
             # Filter the dataframe to include only rows where column 'A' has duplicates
             viral_segments = segmented.loc[duplicates_mask]
             
             # Further filter for only reads containing segments with opposite complementarity
             viral_copybacks = viral_segments.loc[copyback_filter(viral_segments)]
             
    
             # Parse cigar strings to map out reads
             parsed_reads_df, MINDs = melt_segments(viral_copybacks)
             
             plotreads(base_dir, samfilename, reference_name, parsed_reads_df, MINDs)