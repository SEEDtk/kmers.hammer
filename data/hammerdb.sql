CREATE TABLE Hammer (
	hammer VARCHAR(20) PRIMARY KEY,
	genome_id VARCHAR(20) NOT NULL
	);
	INSERT INTO _diagram (table_name, rloc, cloc, description) VALUES ('Hammer', 1, 1, 
		'A hammer is a kmer that is unique to a particular representative group.');
	INSERT INTO _fields (table_name, field_name, field_type, description) VALUES ('Hammer', 'hammer', 'STRING',
		'kmer used to identify a particular representative group');
	INSERT INTO _fields (table_name, field_name, field_type, description) VALUES ('Hammer', 'genome_id', 'STRING',
		'ID of the genome containing the hammer sequence');
