BEGIN TRANSACTION;
-- Table that tracks each pipeline
CREATE TABLE IF NOT EXISTS "pipeline" (
	"pipelineId" integer primary key,
	"statusCode" integer NOT NULL,
	"statusMessage" text,
	"startTimestamp" text,
	"modifiedTimestamp" text,
	"tempFolder" text NOT NULL,
	"outputFolder" text NOT NULL,
	"inputfileBP" integer NOT NULL,
	"info" text
);
-- Table that tracks various scripts run
CREATE TABLE IF NOT EXISTS "scriptUse" (
	"runId"	integer primary key,
	"pipelineId" integer NOT NULL,
	"scriptName" text NOT NULL,
	"start" text integer NOT NULL,
	"end" text,
	"status" text,
	"info" text,
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Save the arguments with which each script runs
CREATE TABLE IF NOT EXISTS "scriptArguments" (
	"argId"	integer primary key,
	"runId"	integer NOT NULL,
	"scriptName" text NOT NULL,
	"argument" text,
	"value" text,
	FOREIGN KEY("runId") REFERENCES "scriptUse"("runId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Table that keeps general logs
CREATE TABLE IF NOT EXISTS "logs" (
    "logId"	integer primary key,
	"runId"	integer NOT NULL,
	"tool" text NOT NULL,
	"timeStamp" integer NOT NULL,
	"actionId" integer NOT NULL,
	"actionName" text,
	FOREIGN KEY("runId") REFERENCES "scriptUse"("runId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Table that stores options used to prepare BLAST
CREATE TABLE IF NOT EXISTS "blastPrepOptions" (
    "prepId" integer primary key,
	"runId"	integer NOT NULL,
	"option" text,
	"value" text,
	FOREIGN KEY("runId") REFERENCES "scriptUse"("runId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Table that stores BLAST submissions and status
CREATE TABLE IF NOT EXISTS "blastSubmissions" (
    "submId" integer primary key,
	"pipelineId"	integer NOT NULL,
	"runId"	integer NOT NULL,
	"RID" text,
	"timeStamp" integer,
	"tempName" text,
	"fastaFile" text,
	"statusCode" integer,
	"statusMessage" text,
	"folder" text,
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId") ON UPDATE CASCADE ON DELETE CASCADE
	FOREIGN KEY("runId") REFERENCES "scriptUse"("runId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Table that stores all antibiotics
CREATE TABLE IF NOT EXISTS "antibiotics" (
    "antibioticId" integer primary key,
	"name" text,
	"class" text
);
-- Table that stores all ARG
CREATE TABLE IF NOT EXISTS "ARG" (
    "geneId" integer primary key,
	"prot" text,
	"nucl" text,
	"gene" test,
	"subtype" text,
	"info" text,
	"nBases" integer,
	"clusterNr" text
);
-- Table that stores all ARG detected in a sample
CREATE TABLE IF NOT EXISTS "detectedARG" (
	"pipelineId" integer NOT NULL,
	"runId" integer NOT NULL,
    "geneId" integer NOT NULL,
	"ARGgroup" integer,
	"fileDepth" real,
	"nSeg" integer,
	"LNsum" integer,
	"KCsum" integer,
	"startPerc" real,
	"LNmax" integer,
	"KCmax" integer,
	"startDepth" real,
	"cover1" real,
	"cover2" real,
	"type" text,
	PRIMARY KEY("runId", "geneId"),
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId") ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY("geneId") REFERENCES "ARG"("geneId") ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY("runId") REFERENCES "scriptUse"("runId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Table that stores all detected bacteria
CREATE TABLE IF NOT EXISTS "detectedBact" (
	"pipelineId" integer NOT NULL,
    "taxId" integer NOT NULL,
	"membership" integer NOT NULL,
	"genus" text,
	"species" text,
	"prob" numeric,
	"relAbundance" numeric,
	"value" numeric,
	PRIMARY KEY("pipelineId", "taxId"),
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Table that stores all AMR predictions
CREATE TABLE IF NOT EXISTS "AMRprediction" (
	"predictionId" integer primary key,
	"pipelineId" integer NOT NULL,
	"membership" integer NOT NULL,
	"antibioticId" integer NOT NULL,
	"resistance" text,
	"value" numeric,
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId") ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY("antibioticId") REFERENCES "antibiotics"("antibioticId") ON UPDATE CASCADE ON DELETE CASCADE
);
-- Add the antibiotics
INSERT INTO antibiotics (name, class)
VALUES
('amikacin', 'Aminoglycosides'),
('amoxicillin', 'Penicillins'),
('amoxicillin-clavulanate', 'Penicillin combinations'),
('ampicillin', 'Penicillins'),
('ampicillin-sulbactam', 'Penicillin combinations'),
('arsphenamine', 'Others'),
('azithromycin', 'Macrolides'),
('azlocillin', 'Penicillins'),
('aztreonam', 'Monobactams'),
('bacitracin', 'Polypeptides'),
('capreomycin', 'Drugs against mycobacteria'),
('cefaclor', 'Cephalosporins (Second generation)'),
('cefadroxil', 'Cephalosporins (First generation)'),
('cefalexin', 'Cephalosporins (First generation)'),
('cefamandole', 'Cephalosporins (Second generation)'),
('cefazolin', 'Cephalosporins (First generation)'),
('cefdinir', 'Cephalosporins (Third generation)'),
('cefditoren', 'Cephalosporins (Third generation)'),
('cefepime', 'Cephalosporins (Fourth generation)'),
('cefixime', 'Cephalosporins (Third generation)'),
('cefmetazole', 'Cephalosporins (Second generation)'),
('cefonicid', 'Cephalosporins (Second generation)'),
('cefoperazone', 'Cephalosporins (Third generation)'),
('cefotaxime', 'Cephalosporins (Third generation)'),
('cefotaxime-clavulanate', 'Penicillin combinations'),
('cefotetan', 'Cephalosporins (Second generation)'),
('cefoxitin', 'Cephalosporins (Second generation)'),
('cefpodoxime', 'Cephalosporins (Third generation)'),
('cefprozil', 'Cephalosporins (Second generation)'),
('ceftaroline fosamil', 'Cephalosporins (Fifth generation)'),
('ceftazidime', 'Cephalosporins (Third generation)'),
('ceftazidime-clavulanate', 'Penicillin combinations'),
('ceftibuten', 'Cephalosporins (Third generation)'),
('ceftizoxime', 'Cephalosporins (Third generation)'),
('ceftobiprole', 'Cephalosporins (Fifth generation)'),
('ceftriaxone', 'Cephalosporins (Third generation)'),
('cefuroxime', 'Cephalosporins (Second generation)'),
('cephalothin', 'Cephalosporins (First generation)'),
('cephapirin', 'Cephalosporins (First generation)'),
('cephradine', 'Cephalosporins (First generation)'),
('chloramphenicol', 'Others'),
('ciprofloxacin', 'Quinolones/Fluoroquinolones'),
('clarithromycin', 'Macrolides'),
('clindamycin', 'Lincosamides'),
('clofazimine', 'Drugs against mycobacteria'),
('colistin', 'Polypeptides'),
('cycloserine', 'Drugs against mycobacteria'),
('dalbavancin', 'Glycopeptides'),
('dapsone', 'Drugs against mycobacteria'),
('daptomycin', 'Lincosamides'),
('demeclocycline', 'Tetracyclines'),
('dicloxacillin', 'Penicillins'),
('doripenem', 'Carbapenems'),
('doxycycline', 'Tetracyclines'),
('enoxacin', 'Quinolones/Fluoroquinolones'),
('ertapenem', 'Carbapenems'),
('erythromycin', 'Macrolides'),
('ethambutol', 'Drugs against mycobacteria'),
('ethionamide', 'Drugs against mycobacteria'),
('flucloxacillin', 'Penicillins'),
('fosfomycin', 'Others'),
('furazolidone', 'Nitrofurans'),
('fusidic acid', 'Others'),
('gatifloxacin', 'Quinolones/Fluoroquinolones'),
('geldanamycin', 'Ansamycins'),
('gemifloxacin', 'Quinolones/Fluoroquinolones'),
('gentamicin', 'Aminoglycosides'),
('grepafloxacin', 'Quinolones/Fluoroquinolones'),
('herbimycin', 'Ansamycins'),
('imipenem', 'Carbapenems'),
('isoniazid', 'Drugs against mycobacteria'),
('kanamycin', 'Aminoglycosides'),
('levofloxacin', 'Quinolones/Fluoroquinolones'),
('lincomycin', 'Lincosamides'),
('linezolid', 'Oxazolidinones'),
('lipopeptide', 'Lincosamides'),
('lomefloxacin', 'Quinolones/Fluoroquinolones'),
('loracarbef', 'Carbacephem'),
('mafenide', 'Sulfonamides'),
('meropenem', 'Carbapenems'),
('metacycline', 'Tetracyclines'),
('methicillin', 'Penicillins'),
('metronidazole', 'Others'),
('mezlocillin', 'Penicillins'),
('minocycline', 'Tetracyclines'),
('moxalactam', 'Cephalosporins (Third generation)'),
('moxifloxacin', 'Quinolones/Fluoroquinolones'),
('mupirocin', 'Others'),
('nadifloxacin', 'Quinolones/Fluoroquinolones'),
('nafcillin', 'Penicillins'),
('nalidixic acid', 'Quinolones/Fluoroquinolones'),
('neomycin', 'Aminoglycosides'),
('netilmicin', 'Aminoglycosides'),
('nitrofurantoin', 'Nitrofurans'),
('norfloxacin', 'Quinolones/Fluoroquinolones'),
('ofloxacin', 'Quinolones/Fluoroquinolones'),
('oritavancin', 'Glycopeptides'),
('oxacillin', 'Penicillins'),
('oxytetracycline', 'Tetracyclines'),
('paromomycin', 'Aminoglycosides'),
('penicillin g', 'Penicillins'),
('penicillin v', 'Penicillins'),
('piperacillin', 'Penicillins'),
('piperacillin-tazobactam', 'Penicillin combinations'),
('platensimycin', 'Others'),
('polymyxin b', 'Polypeptides'),
('posizolid', 'Oxazolidinones'),
('pyrazinamide', 'Drugs against mycobacteria'),
('quinupristin-dalfopristin', 'Others'),
('radezolid', 'Oxazolidinones'),
('rifabutin', 'Drugs against mycobacteria'),
('rifampicin', 'Drugs against mycobacteria'),
('rifapentine', 'Drugs against mycobacteria'),
('rifaximin', 'Ansamycins'),
('roxithromycin', 'Macrolides'),
('silver sulfadiazine', 'Sulfonamides'),
('sparfloxacin', 'Quinolones/Fluoroquinolones'),
('spectinomycin', 'Aminoglycosides'),
('spiramycin', 'Macrolides'),
('streptomycin', 'Aminoglycosides'),
('sulfacetamide', 'Sulfonamides'),
('sulfadiazine', 'Sulfonamides'),
('sulfadimethoxine', 'Sulfonamides'),
('sulfamethizole', 'Sulfonamides'),
('sulfamethoxazole', 'Sulfonamides'),
('sulfasalazine', 'Sulfonamides'),
('sulfisoxazole', 'Sulfonamides'),
('teicoplanin', 'Glycopeptides'),
('telavancin', 'Glycopeptides'),
('telithromycin', 'Macrolides'),
('temafloxacin', 'Quinolones/Fluoroquinolones'),
('temocillin', 'Penicillins'),
('tetracycline', 'Tetracyclines'),
('thiamphenicol', 'Others'),
('ticarcillin', 'Penicillins'),
('ticarcillin-clavulanate', 'Penicillin combinations'),
('tigecycline', 'Others'),
('tinidazole', 'Others'),
('tobramycin', 'Aminoglycosides'),
('torezolid', 'Oxazolidinones'),
('trimethoprim', 'Others'),
('trimethoprim-sulfamethoxazole', 'Sulfonamides'),
('trovafloxacin', 'Quinolones/Fluoroquinolones'),
('unknown', 'Unknown'),
('vancomycin', 'Glycopeptides');

COMMIT;
