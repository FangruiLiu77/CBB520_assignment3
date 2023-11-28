import os
import csv
import subprocess
from flask import Flask, request, render_template, redirect, url_for
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from werkzeug.utils import secure_filename


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = '/home/fl118/CBB520_assignment3/upload'

# Path to your database.fasta file
database_path = "/home/fl118/CBB520_assignment3/database"

# Load the database sequences
#database_sequences = list(SeqIO.parse(database_path, "fasta"))

@app.route('/', methods=['GET'])
def index():
    
    return render_template('my_web.html')

@app.route('/blast-results', methods=['POST'])
def submit():
    user_input = request.form['user_input']

    # Create a temporary query sequence file
    query_seq_file = "user_query.fasta"
    with open(query_seq_file, "w") as file:
        file.write(f">User_Query\n{user_input}")

    summary = []
    user_input_pI = calculate_pI(user_input)
    # isoelectric_point.append({"User Input Isoelectric Point": f"{user_input_pI:.2f}"})
    with open("database.fasta", 'r') as data:
        lines = data.readlines()
        idx = [lines[i] for i in range(0, len(lines), 2)]
        iso = [calculate_pI(lines[i]) for i in range(1, len(lines)+1, 2)]
    
    isoelectrics = [{idx: iso} for idx, iso in zip(idx, iso)]
            

    # Run a BLAST search using BLAST+ and parse the result
    out_filename = "blast_results.xml"
    command = f'blastp -query {query_seq_file} -db {database_path} -out {out_filename} -outfmt 5 -max_target_seqs 10'
    subprocess.run(command, shell=True)
    with open(out_filename, "r") as blast_output:
        blast_records = NCBIXML.parse(blast_output)
        for record in blast_records:
            for alignment in record.alignments:
                identity = (alignment.hsps[0].identities / alignment.hsps[0].align_length) * 100
                similarity = (alignment.hsps[0].positives / alignment.hsps[0].align_length) * 100

                summary.append({
                    "Hit ID": alignment.hit_id,
                    "Identity": f"{identity:.2f}% identical",
                    "Similarity": f"{similarity:.2f}% similar"
                })

    os.remove(query_seq_file)
    os.remove(out_filename)

    #return jsonify({"summary": jsonify(summary), "isoelectrics": jsonify(isoelectrics)})
    return render_template('blast_results.html', isoelectrics=isoelectrics, 
                           user_input_pI=user_input_pI, 
                           summary=summary, )


def calculate_pI(sequence):
    # Calculate the isoelectric point using ProtParam
    protein_analyzer = ProtParam.ProteinAnalysis(sequence)
    return protein_analyzer.isoelectric_point()


snp_data = {}

# Specify the path to your CSV file
csv_file = '/home/fl118/CBB520_assignment3/info.csv'

# Open and read the CSV file
with open(csv_file, 'r') as file:
    csv_reader = csv.DictReader(file)
    for row in csv_reader:
        dsSNP_number = row['dsSNP_number']
        snp_data[dsSNP_number] = {
            'dsSNP_number': dsSNP_number,
            'intron': row['intron'],
            'coding_region': row['coding_region'],
            'synonymous':  row['synonymous '],
            'utr5': row['utr5'],
            'utr3': row['utr3'],
            'upstream': row['upstream'],
            'downstream': row['downstream'],
            'Ref_Allele': row['Ref_Allele'],
            'Alt_Allele': row['Alt_Allele'],
            'dbsnp_link': row['dbsnp_link'],
            'omim_link': row['omim_link']
        }

@app.route('/snp_by_number', methods=['POST'])
def snp_by_number():
    snp_number = request.form['user_snp_number']
    if snp_number in snp_data:
        result = snp_data[snp_number]
    else:
        result = 'This dsSNP number is not in our database.'
   
    return render_template('snp_by_number.html', result=result)


@app.route('/snp_show_all', methods=['GET'])
def snp_show_all():
    return render_template('snp_show_all.html', snp_data= snp_data)


@app.route('/snp_show_synonymous', methods=['GET'])
def snp_show_synonymous():
    return render_template('snp_show_synonymous.html', snp_data= snp_data)


@app.route('/snp_show_non_synonymous', methods=['GET'])
def snp_show_non_synonymous():
    return render_template('snp_show_non_synonymous.html', snp_data= snp_data)


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/upload_pdb', methods=['GET', 'POST'])
def upload_pdb():
    if request.method == 'POST':
        # Check if the post request has the file part
        if 'pdb_file' not in request.files:
            return redirect(request.url)
        file = request.files['pdb_file']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            return redirect(request.url)
        if file:
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return redirect(url_for('index'))
    return render_template('upload_pdb.html')


def run_tmalign(pdb_file1, pdb_file2):
    command = f"TMalign {pdb_file1} {pdb_file2}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    output, error = process.communicate()

    if error:
        print("Error in running TMalign: ", error)
    return output.decode()

@app.route('/run_alignment', methods=['POST'])
def run_alignment():
    pdb_file1 = request.form.get('pdb1')
    pdb_file2 = request.form.get('pdb2')
    
    alignment_result = run_tmalign(pdb_file1, pdb_file2)
    # Process and return the alignment result
    return render_template('alignment_result.html', result=alignment_result)

if __name__ == '__main__':
    app.run(host="0.0.0.0", port="49151")
