#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import sys
import subprocess
import pandas as pd
import multiprocessing
from Bio import AlignIO
from pprint import pprint
from functools import partial

def getListTree():
    listTree = ['grpOrtho_10001',  'grpOrtho_10013',  'grpOrtho_10016',  'grpOrtho_10022',  'grpOrtho_10027',  'grpOrtho_10030',  'grpOrtho_10031',  'grpOrtho_10033',  'grpOrtho_10034',  'grpOrtho_10039',  'grpOrtho_10040',  'grpOrtho_10041',  'grpOrtho_10042',  'grpOrtho_10050',  'grpOrtho_10051',  'grpOrtho_10052',  'grpOrtho_10056',  'grpOrtho_10058',  'grpOrtho_10059',  'grpOrtho_10066',  'grpOrtho_10068',  'grpOrtho_10069',  'grpOrtho_10076',  'grpOrtho_10089',  'grpOrtho_100',  'grpOrtho_10102',  'grpOrtho_10113',  'grpOrtho_10114',  'grpOrtho_10118',  'grpOrtho_10124',  'grpOrtho_10128',  'grpOrtho_10129',  'grpOrtho_10130',  'grpOrtho_10136',  'grpOrtho_10141',  'grpOrtho_10143',  'grpOrtho_10144',  'grpOrtho_10147',  'grpOrtho_10149',  'grpOrtho_10150',  'grpOrtho_10153',  'grpOrtho_10156',  'grpOrtho_10159',  'grpOrtho_10161',  'grpOrtho_10166',  'grpOrtho_10167',  'grpOrtho_10169',  'grpOrtho_10172',  'grpOrtho_10174',  'grpOrtho_10176',  'grpOrtho_10177',  'grpOrtho_10182',  'grpOrtho_10192',  'grpOrtho_10194',  'grpOrtho_10196',  'grpOrtho_10197',  'grpOrtho_10198',  'grpOrtho_10215',  'grpOrtho_10217',  'grpOrtho_10223',  'grpOrtho_10226',  'grpOrtho_10227',  'grpOrtho_10230',  'grpOrtho_10231',  'grpOrtho_10232',  'grpOrtho_1023',  'grpOrtho_10340',  'grpOrtho_10383',  'grpOrtho_10403',  'grpOrtho_10634',  'grpOrtho_1063',  'grpOrtho_10',  'grpOrtho_11076',  'grpOrtho_111',  'grpOrtho_113',  'grpOrtho_11',  'grpOrtho_1273',  'grpOrtho_13111',  'grpOrtho_13114',  'grpOrtho_13115',  'grpOrtho_13118',  'grpOrtho_13119',  'grpOrtho_13120',  'grpOrtho_13129',  'grpOrtho_13134',  'grpOrtho_13139',  'grpOrtho_13142',  'grpOrtho_13144',  'grpOrtho_13145',  'grpOrtho_13151',  'grpOrtho_13152',  'grpOrtho_13159',  'grpOrtho_13162',  'grpOrtho_13167',  'grpOrtho_13171',  'grpOrtho_13174',  'grpOrtho_13175',  'grpOrtho_13176',  'grpOrtho_13179',  'grpOrtho_13180',  'grpOrtho_13182',  'grpOrtho_13183',  'grpOrtho_13184',  'grpOrtho_13187',  'grpOrtho_13192',  'grpOrtho_13193',  'grpOrtho_13194',  'grpOrtho_13208',  'grpOrtho_13216',  'grpOrtho_13223',  'grpOrtho_13225',  'grpOrtho_13229',  'grpOrtho_13232',  'grpOrtho_13233',  'grpOrtho_13235',  'grpOrtho_13236',  'grpOrtho_13237',  'grpOrtho_13244',  'grpOrtho_13245',  'grpOrtho_13246',  'grpOrtho_13248',  'grpOrtho_13251',  'grpOrtho_13257',  'grpOrtho_13258',  'grpOrtho_13260',  'grpOrtho_13261',  'grpOrtho_13262',  'grpOrtho_13267',  'grpOrtho_13269',  'grpOrtho_13271',  'grpOrtho_13273',  'grpOrtho_13275',  'grpOrtho_13282',  'grpOrtho_13287',  'grpOrtho_13292',  'grpOrtho_13293',  'grpOrtho_13297',  'grpOrtho_13298',  'grpOrtho_13300',  'grpOrtho_13306',  'grpOrtho_13311',  'grpOrtho_13313',  'grpOrtho_13315',  'grpOrtho_13322',  'grpOrtho_13323',  'grpOrtho_13325',  'grpOrtho_13326',  'grpOrtho_13332',  'grpOrtho_13340',  'grpOrtho_13341',  'grpOrtho_13343',  'grpOrtho_13351',  'grpOrtho_13354',  'grpOrtho_13355',  'grpOrtho_13356',  'grpOrtho_13359',  'grpOrtho_13367',  'grpOrtho_13368',  'grpOrtho_13371',  'grpOrtho_13384',  'grpOrtho_13385',  'grpOrtho_13387',  'grpOrtho_13391',  'grpOrtho_13392',  'grpOrtho_13403',  'grpOrtho_13404',  'grpOrtho_13407',  'grpOrtho_13419',  'grpOrtho_13420',  'grpOrtho_13422',  'grpOrtho_13423',  'grpOrtho_13426',  'grpOrtho_13429',  'grpOrtho_13430',  'grpOrtho_13433',  'grpOrtho_13438',  'grpOrtho_13439',  'grpOrtho_13441',  'grpOrtho_13447',  'grpOrtho_13458',  'grpOrtho_13472',  'grpOrtho_13473',  'grpOrtho_13476',  'grpOrtho_13480',  'grpOrtho_13485',  'grpOrtho_13487',  'grpOrtho_13495',  'grpOrtho_13498',  'grpOrtho_13501',  'grpOrtho_13502',  'grpOrtho_13516',  'grpOrtho_13521',  'grpOrtho_13522',  'grpOrtho_13526',  'grpOrtho_13529',  'grpOrtho_13530',  'grpOrtho_13541',  'grpOrtho_13542',  'grpOrtho_13543',  'grpOrtho_13546',  'grpOrtho_13547',  'grpOrtho_13630',  'grpOrtho_13645',  'grpOrtho_136',  'grpOrtho_13799',  'grpOrtho_13823',  'grpOrtho_13824',  'grpOrtho_13828',  'grpOrtho_13829',  'grpOrtho_13838',  'grpOrtho_13877',  'grpOrtho_13885',  'grpOrtho_13886',  'grpOrtho_13887',  'grpOrtho_13893',  'grpOrtho_13906',  'grpOrtho_13910',  'grpOrtho_13925',  'grpOrtho_13929',  'grpOrtho_13933',  'grpOrtho_13936',  'grpOrtho_13944',  'grpOrtho_13949',  'grpOrtho_13950',  'grpOrtho_13952',  'grpOrtho_13966',  'grpOrtho_13976',  'grpOrtho_13981',  'grpOrtho_139',  'grpOrtho_13',  'grpOrtho_14012',  'grpOrtho_14030',  'grpOrtho_14035',  'grpOrtho_14050',  'grpOrtho_14057',  'grpOrtho_14060',  'grpOrtho_14071',  'grpOrtho_14082',  'grpOrtho_14085',  'grpOrtho_140',  'grpOrtho_14109',  'grpOrtho_14119',  'grpOrtho_14123',  'grpOrtho_14136',  'grpOrtho_14153',  'grpOrtho_14162',  'grpOrtho_14169',  'grpOrtho_14173',  'grpOrtho_14177',  'grpOrtho_14187',  'grpOrtho_14191',  'grpOrtho_14193',  'grpOrtho_14214',  'grpOrtho_14217',  'grpOrtho_14219',  'grpOrtho_14220',  'grpOrtho_14223',  'grpOrtho_14226',  'grpOrtho_14231',  'grpOrtho_14245',  'grpOrtho_14248',  'grpOrtho_14255',  'grpOrtho_14259',  'grpOrtho_14260',  'grpOrtho_14281',  'grpOrtho_14297',  'grpOrtho_14313',  'grpOrtho_14324',  'grpOrtho_14346',  'grpOrtho_14362',  'grpOrtho_14365',  'grpOrtho_14383',  'grpOrtho_14395',  'grpOrtho_14397',  'grpOrtho_14446',  'grpOrtho_14448',  'grpOrtho_14449',  'grpOrtho_14456',  'grpOrtho_14457',  'grpOrtho_14484',  'grpOrtho_14487',  'grpOrtho_14492',  'grpOrtho_14501',  'grpOrtho_14520',  'grpOrtho_14534',  'grpOrtho_14536',  'grpOrtho_14544',  'grpOrtho_1454',  'grpOrtho_14590',  'grpOrtho_14599',  'grpOrtho_145',  'grpOrtho_14603',  'grpOrtho_14618',  'grpOrtho_14643',  'grpOrtho_14662',  'grpOrtho_14673',  'grpOrtho_14699',  'grpOrtho_14762',  'grpOrtho_14798',  'grpOrtho_147',  'grpOrtho_14800',  'grpOrtho_14829',  'grpOrtho_14985',  'grpOrtho_15138',  'grpOrtho_15270',  'grpOrtho_15330',  'grpOrtho_15341',  'grpOrtho_15378',  'grpOrtho_1542',  'grpOrtho_15496',  'grpOrtho_15525',  'grpOrtho_15558',  'grpOrtho_15630',  'grpOrtho_15654',  'grpOrtho_15699',  'grpOrtho_15708',  'grpOrtho_15709',  'grpOrtho_15710',  'grpOrtho_15714',  'grpOrtho_15715',  'grpOrtho_15726',  'grpOrtho_1572',  'grpOrtho_15758',  'grpOrtho_15775',  'grpOrtho_15782',  'grpOrtho_15784',  'grpOrtho_15785',  'grpOrtho_15845',  'grpOrtho_15878',  'grpOrtho_15946',  'grpOrtho_1596',  'grpOrtho_15984',  'grpOrtho_16292',  'grpOrtho_16343',  'grpOrtho_16535',  'grpOrtho_16631',  'grpOrtho_16731',  'grpOrtho_16772',  'grpOrtho_167',  'grpOrtho_16831',  'grpOrtho_1705',  'grpOrtho_171',  'grpOrtho_17206',  'grpOrtho_1720',  'grpOrtho_17536',  'grpOrtho_1792',  'grpOrtho_17932',  'grpOrtho_1805',  'grpOrtho_18229',  'grpOrtho_1846',  'grpOrtho_18476',  'grpOrtho_18849',  'grpOrtho_1900',  'grpOrtho_191',  'grpOrtho_19311',  'grpOrtho_1935',  'grpOrtho_19494',  'grpOrtho_19644',  'grpOrtho_196',  'grpOrtho_1970',  'grpOrtho_198',  'grpOrtho_20015',  'grpOrtho_20070',  'grpOrtho_202',  'grpOrtho_20328',  'grpOrtho_203',  'grpOrtho_2041',  'grpOrtho_20533',  'grpOrtho_205',  'grpOrtho_2060',  'grpOrtho_206',  'grpOrtho_2070',  'grpOrtho_20816',  'grpOrtho_20954',  'grpOrtho_20956',  'grpOrtho_20990',  'grpOrtho_21178',  'grpOrtho_21271',  'grpOrtho_212',  'grpOrtho_2139',  'grpOrtho_214',  'grpOrtho_2170',  'grpOrtho_220',  'grpOrtho_225',  'grpOrtho_227',  'grpOrtho_228',  'grpOrtho_22',  'grpOrtho_230',  'grpOrtho_2314',  'grpOrtho_23',  'grpOrtho_2401',  'grpOrtho_244',  'grpOrtho_2464',  'grpOrtho_2481',  'grpOrtho_2532',  'grpOrtho_2533',  'grpOrtho_2548',  'grpOrtho_2606',  'grpOrtho_2613',  'grpOrtho_2626',  'grpOrtho_264',  'grpOrtho_266',  'grpOrtho_26714',  'grpOrtho_2767',  'grpOrtho_27',  'grpOrtho_281',  'grpOrtho_284',  'grpOrtho_286',  'grpOrtho_2885',  'grpOrtho_288',  'grpOrtho_28',  'grpOrtho_293',  'grpOrtho_29',  'grpOrtho_300',  'grpOrtho_3018',  'grpOrtho_301',  'grpOrtho_3020',  'grpOrtho_302',  'grpOrtho_3046',  'grpOrtho_3056',  'grpOrtho_305',  'grpOrtho_3064',  'grpOrtho_307',  'grpOrtho_311',  'grpOrtho_312',  'grpOrtho_314',  'grpOrtho_329',  'grpOrtho_332',  'grpOrtho_339',  'grpOrtho_342',  'grpOrtho_344',  'grpOrtho_3451',  'grpOrtho_352',  'grpOrtho_354',  'grpOrtho_357',  'grpOrtho_362',  'grpOrtho_3662',  'grpOrtho_369',  'grpOrtho_3729',  'grpOrtho_372',  'grpOrtho_3738',  'grpOrtho_3748',  'grpOrtho_3765',  'grpOrtho_3768',  'grpOrtho_379',  'grpOrtho_380',  'grpOrtho_38',  'grpOrtho_3900',  'grpOrtho_3903',  'grpOrtho_3909',  'grpOrtho_3913',  'grpOrtho_392',  'grpOrtho_398',  'grpOrtho_4054',  'grpOrtho_408',  'grpOrtho_411',  'grpOrtho_413',  'grpOrtho_4191',  'grpOrtho_4195',  'grpOrtho_419',  'grpOrtho_424',  'grpOrtho_4272',  'grpOrtho_427',  'grpOrtho_4280',  'grpOrtho_4341',  'grpOrtho_4364',  'grpOrtho_4367',  'grpOrtho_438',  'grpOrtho_439',  'grpOrtho_4415',  'grpOrtho_4432',  'grpOrtho_4479',  'grpOrtho_453',  'grpOrtho_4542',  'grpOrtho_456',  'grpOrtho_4570',  'grpOrtho_458',  'grpOrtho_461',  'grpOrtho_4677',  'grpOrtho_469',  'grpOrtho_4723',  'grpOrtho_482',  'grpOrtho_484',  'grpOrtho_4867',  'grpOrtho_486',  'grpOrtho_489',  'grpOrtho_4927',  'grpOrtho_494',  'grpOrtho_4964',  'grpOrtho_497',  'grpOrtho_501',  'grpOrtho_502',  'grpOrtho_5093',  'grpOrtho_5131',  'grpOrtho_5183',  'grpOrtho_5234',  'grpOrtho_5236',  'grpOrtho_5242',  'grpOrtho_5254',  'grpOrtho_5328',  'grpOrtho_5481',  'grpOrtho_5563',  'grpOrtho_5601',  'grpOrtho_560',  'grpOrtho_5658',  'grpOrtho_56',  'grpOrtho_5719',  'grpOrtho_5723',  'grpOrtho_5778',  'grpOrtho_5814',  'grpOrtho_5842',  'grpOrtho_5863',  'grpOrtho_5879',  'grpOrtho_5903',  'grpOrtho_5980',  'grpOrtho_6110',  'grpOrtho_6216',  'grpOrtho_6270',  'grpOrtho_62',  'grpOrtho_6309',  'grpOrtho_6475',  'grpOrtho_6497',  'grpOrtho_6498',  'grpOrtho_680',  'grpOrtho_6816',  'grpOrtho_6844',  'grpOrtho_6856',  'grpOrtho_6860',  'grpOrtho_6876',  'grpOrtho_68',  'grpOrtho_6953',  'grpOrtho_7004',  'grpOrtho_7021',  'grpOrtho_7042',  'grpOrtho_7115',  'grpOrtho_7142',  'grpOrtho_7145',  'grpOrtho_7187',  'grpOrtho_7229',  'grpOrtho_7250',  'grpOrtho_7261',  'grpOrtho_72',  'grpOrtho_7389',  'grpOrtho_746',  'grpOrtho_7488',  'grpOrtho_74',  'grpOrtho_7514',  'grpOrtho_7691',  'grpOrtho_7853',  'grpOrtho_7861',  'grpOrtho_7875',  'grpOrtho_7879',  'grpOrtho_7882',  'grpOrtho_7883',  'grpOrtho_7885',  'grpOrtho_7886',  'grpOrtho_7889',  'grpOrtho_7899',  'grpOrtho_7902',  'grpOrtho_7911',  'grpOrtho_7916',  'grpOrtho_7917',  'grpOrtho_7940',  'grpOrtho_7952',  'grpOrtho_7955',  'grpOrtho_7957',  'grpOrtho_7962',  'grpOrtho_7979',  'grpOrtho_7982',  'grpOrtho_7984',  'grpOrtho_7987',  'grpOrtho_7991',  'grpOrtho_8005',  'grpOrtho_8011',  'grpOrtho_8013',  'grpOrtho_8016',  'grpOrtho_8040',  'grpOrtho_8053',  'grpOrtho_8099',  'grpOrtho_8106',  'grpOrtho_8107',  'grpOrtho_8110',  'grpOrtho_8135',  'grpOrtho_8136',  'grpOrtho_8202',  'grpOrtho_8207',  'grpOrtho_8241',  'grpOrtho_8390',  'grpOrtho_8483',  'grpOrtho_84',  'grpOrtho_8503',  'grpOrtho_8514',  'grpOrtho_8534',  'grpOrtho_8626',  'grpOrtho_8629',  'grpOrtho_8645',  'grpOrtho_8646',  'grpOrtho_8658',  'grpOrtho_8697',  'grpOrtho_8701',  'grpOrtho_8703',  'grpOrtho_8711',  'grpOrtho_8798',  'grpOrtho_8831',  'grpOrtho_8845',  'grpOrtho_8861',  'grpOrtho_88',  'grpOrtho_89',  'grpOrtho_9032',  'grpOrtho_9043',  'grpOrtho_9130',  'grpOrtho_9161',  'grpOrtho_9195',  'grpOrtho_9214',  'grpOrtho_9248',  'grpOrtho_9255',  'grpOrtho_9270',  'grpOrtho_9293',  'grpOrtho_9301',  'grpOrtho_9307',  'grpOrtho_9313',  'grpOrtho_9325',  'grpOrtho_940',  'grpOrtho_9431',  'grpOrtho_9449',  'grpOrtho_9458',  'grpOrtho_9461',  'grpOrtho_9466',  'grpOrtho_94',  'grpOrtho_9500',  'grpOrtho_9529',  'grpOrtho_9551',  'grpOrtho_9555',  'grpOrtho_9566',  'grpOrtho_9593',  'grpOrtho_9645',  'grpOrtho_9660',  'grpOrtho_9682',  'grpOrtho_9699',  'grpOrtho_9767',  'grpOrtho_9812',  'grpOrtho_9815',  'grpOrtho_9817',  'grpOrtho_9819',  'grpOrtho_9828',  'grpOrtho_9829',  'grpOrtho_9833',  'grpOrtho_9834',  'grpOrtho_9844',  'grpOrtho_9851',  'grpOrtho_9857',  'grpOrtho_9867',  'grpOrtho_9869',  'grpOrtho_9874',  'grpOrtho_9879',  'grpOrtho_9881',  'grpOrtho_9882',  'grpOrtho_9884',  'grpOrtho_9886',  'grpOrtho_9887',  'grpOrtho_9888',  'grpOrtho_9893',  'grpOrtho_9895',  'grpOrtho_9897',  'grpOrtho_9898',  'grpOrtho_9900',  'grpOrtho_9904',  'grpOrtho_9905',  'grpOrtho_9913',  'grpOrtho_9914',  'grpOrtho_9915',  'grpOrtho_9919',  'grpOrtho_9927',  'grpOrtho_9929',  'grpOrtho_9936',  'grpOrtho_9939',  'grpOrtho_9949',  'grpOrtho_9950',  'grpOrtho_9951',  'grpOrtho_9953',  'grpOrtho_9955',  'grpOrtho_9958',  'grpOrtho_9964',  'grpOrtho_9969',  'grpOrtho_9976',  'grpOrtho_9977',  'grpOrtho_9980',  'grpOrtho_9982',  'grpOrtho_9988',  'grpOrtho_9996',  'grpOrtho_9']
    return listTree

def main(path, t):

    res = {'idRegion': [], 'G4NN': [], 'G4H': [], 'cGcC': [], 'qgrs': [], 'gene': []}
    try:
        alignment = AlignIO.read(open(path+'AlignmentOrtho/'+t+'.fa'), 'fasta')
        # alignment = AlignIO.read(open(path+'AlignmentOrtho/Euk/Alignment_11076.fa'), 'fasta')
    except:
        print('No alignment for this tree')
    else:
        with open(path+'mergedpG4FamiliesPecentrOrtho/'+tree+'.csv') as f:
            content = f.read()
            for l in content.split('\n'):
                l = l.split('\t')
                if l[0] != 'Family' and l[0] != '':
                    start = int(l[1])
                    end = int(l[2])

                    for record in alignment:
                        res['idRegion'].append(tree+':'+str(start)+'-'+str(end))
                        res['gene'].append(record.id)
                        res['qgrs'].append(0)
                        res['G4NN'].append(0)
                        res['G4H'].append(0)
                        res['cGcC'].append(0)

                        if str(record.seq[start:end]).translate(None, '-') != '' and len(str(record.seq[start:end]).translate(None, '-')) >= 20:
                            if tree+':'+str(start)+'-'+str(end) == '424:2679-2946':
                                print('>'+record.id)
                                print(str(record.seq[start-1:end]).translate(None, '-'))
                            output = open(path+'tmp.fa', "w")
                            output.write('>'+record.id+'\n')
                            output.write(str(record.seq[start-1:end]).translate(None, '-')+'\n')
                            output.close()

                            qrgs = subprocess.check_output("/home/vana2406/software/qgrs-cpp/qgrs -i "+path+"tmp.fa", shell=True)
                            # qrgs = subprocess.check_output("/home/vana2406/software/qgrs-cpp/qgrs -i "+path+"tmp.fa", shell=True)
                            qrgs = qrgs.decode("utf-8")
                            qrgsRes = qrgs.split('\n')
                            try:
                                qrgsRes.remove('ID          T1  T2  T3  T4 TS  GS  SEQ')
                            except:
                                try:
                                    qrgsRes.remove('ID         T1 T2 T3 T4 TS  GS  SEQ')
                                except:
                                    try:
                                        qrgsRes.remove('ID        T1T2T3T4 TS  GS  SEQ')
                                    except:
                                        try:
                                            qrgsRes.remove('ID           T1   T2   T3   T4 TS  GS  SEQ')
                                        except:
                                            print(qrgsRes)

                            qrgsRes.remove('')
                            qrgsRes.remove('--------------------------------------------------------------------------------------------')

                            if len(qrgsRes)-1 > 0:
                                res['qgrs'][-1] = 1

                            g4rnascreener = subprocess.check_output("/home/vana2406/software/"+\
                            "g4rna_screener/screen.py "+path+"tmp.fa -a /home/vana2406/"+\
                            "software/g4rna_screener/G4RNA_2016-11-07.pkl -w 60 -s 10 "+\
                            "-c description cGcC G4H G4NN sequence start end -e", shell=True)
                            g4rnascreener = g4rnascreener.decode("utf-8")
                            g4rnascreener = g4rnascreener.split('\n')
                            for w in g4rnascreener:
                                w = w.split('\t')
                                if w[0] != '':
                                    if 'FBgn0003345' in w[1] and start == 36725 and end == 36988:
                                        print(w)
                                    # print(w)
                                    if float(w[2]) >= 4.5:
                                        res['cGcC'][-1] = 1
                                    if float(w[3]) >= 0.9:
                                        res['G4H'][-1] = 1
                                    if float(w[7]) >= 0.5:
                                        res['G4NN'][-1] = 1
        tmp = pd.DataFrame(data=res).drop_duplicates(subset=None, keep='first', inplace=False)
        tmp.to_csv(path_or_buf=path+'FamG4Pred/'+tree+'.csv', header=True, index=None, sep='\t')

if __name__ == '__main__':
    path = sys.argv[1]
    listTree = getListTree()

    p = multiprocessing.Pool(15)
	func = partial(main, path)
	data = p.map(func, [ i for i in listTree ])
	p.close
