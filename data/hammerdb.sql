CREATE TABLE Hammer (
    hammer VARCHAR(20) PRIMARY KEY,
    fid VARCHAR(30) NOT NULL,
    strength DOUBLE NOT NULL
    );
    CREATE INDEX idxFid ON Hammer (fid);
    INSERT INTO _diagram (table_name, rloc, cloc, description) VALUES ('Hammer', 1, 1,
        'A hammer is a kmer that is unique to a particular representative group.');
    INSERT INTO _fields (table_name, field_name, field_type, description) VALUES ('Hammer', 'hammer', 'STRING',
        'kmer used to identify a particular representative group');
    INSERT INTO _fields (table_name, field_name, field_type, description) VALUES ('Hammer', 'fid', 'STRING',
        'ID of the feature in the genome containing the hammer sequence');
    INSERT INTO _fields (table_name, field_name, field_type, description) VALUES ('Hammer', 'strength', 'DOUBLE',
        'evidentiary strength of the hammer (ranking for 0 to 1, with 1 indicating best');

