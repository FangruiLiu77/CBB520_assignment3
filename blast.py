import os
from flask import Flask, request, render_template, jsonify
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import ProtParam

app = Flask(__name__)

@app.route('/', methods=['GET'])
def home():
    return render_template("my_web.html")


# Path to your database.fasta file
database_path = "database.fasta"

# Load the database sequences
database_sequences = list(SeqIO.parse(database_path, "fasta"))


@app.route('/submit', methods=['POST'])
def submit():
    user_input = request.form['user_input']
    return user_input

if __name__ == '__main__':
    app.run(port="8080")