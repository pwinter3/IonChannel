import csv, os

def get_seqs():
    upid_list = []
    ionf = open('betse_channel_assignments/data/go/ion-channels.csv')
    ionr = csv.reader(ionf)
    ionr.next()
    for row in ionr:
        upid = row[4]
        upid_list.append(upid)
    ionf.close()
    return upid_list

def download_seqs(upid_list):
    os.chdir('betse_channel_assignments/data/uniprot_seq')
    for upid in upid_list:
        os.system('wget http://www.uniprot.org/uniprot/%s.fasta' % upid)
    os.system('cat *.fasta >../uniprot_seq.fasta')

def makeblastdb():
    os.chdir('betse_channel_assignments/data/uniprot_betse_seq')
    os.syste('formatdb -i uniprot_betse_seq.fasta -p T -o F')
    
def do_blast(upid_list):
    for upid in upid_list[0:1]:
        query_fasta = 'betse_channel_assignments/data/uniprot_seq/%s.fasta' % upid
        os.system('blastall -p blastp ' \
                '-d betse_channel_assignments/data/uniprot_betse_seq/uniprot_betse_seq.fasta' \
                ' -i betse_channel_assignments/data/uniprot_seq/%s.fasta -o betse_channel_assignments/data/blast_res/out.%s' %  \
                (upid, upid))

def do_blastall():
    query_fasta = ''
    os.system('blastall -p blastp -m 9 -e 1e-9' \
            ' -d betse_channel_assignments/data/uniprot_betse_seq/uniprot_betse_seq.fasta' \
            ' -ibetse_channel_assignments/data/uniprot_seq.fasta -o betse_channel_assignments/data/out.blast_res')

def write_matches():
    fbetse = open('betse_channel_assignments/data/betse_channels.csv', 'rU')
    rbetse = csv.reader(fbetse)
    rbetse.next()
    betse_class_lookup = {}
    for row in rbetse:
        betse_module = row[0].split('/')[-1].split('.')[0]
        betse_class = row[1]
        uacc = row[5]
        betse_class_lookup[uacc] = betse_module + '.' + betse_class
    fbetse.close()
    data_dict = {}
    fin = open('betse_channel_assignments/data/out.blast_res')
    for line in fin:
        if line.startswith('#'):
            continue
        line_split = line.split('\t')
        query = line_split[0]
        query = query.split('|')[1]
        target = line_split[1]
        target = target.split('|')[1]
        if not query in data_dict:
            data_dict[query] = []
            data_dict[query].append(target)
        else:
            pass
    fin.close()
    fout = open('betse_channel_assignments/data/matches.csv', 'w')
    for query in data_dict:
        fout.write(query)
        for target in data_dict[query]:
            fout.write(',')
            fout.write(target)
        fout.write('\n')
    fout.close()
    fchannels = open('betse_channel_assignments/data/go/ion-channels.csv')
    rchannels = csv.reader(fchannels)
    row = rchannels.next()
    ffinal = open('betse_channel_assignments/data/assigned_channels/assigned-ion-channels.csv', 'w')
    wfinal = csv.writer(ffinal)
    wfinal.writerow(row + ['BETSE Model'])
    for row in rchannels:
        upid = row[4]
        if upid in data_dict:
            model_code_list = []
            for model in data_dict[upid]:
                model_code = betse_class_lookup[model]
                model_code_list.append(model_code)
            model_str = '|'.join(model_code_list)
            if row[7].startswith('NOT AN ION CHANNEL'):
                print 'matched a non-ion channel'
        else:
            model_str = '-'
        wfinal.writerow(row + [model_str])
    ffinal.close()
    fchannels.close()

get_seqs()
download_seqs()
makeblastdb()
do_blast()
do_blastall()
write_matches()
