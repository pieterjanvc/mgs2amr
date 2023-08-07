BEGIN TRANSACTION;
-- Table 
CREATE TABLE IF NOT EXISTS "pipeline" (
	"pipelineId" integer PRIMARY KEY AUTOINCREMENT,
	"statusCode" integer,
	"name" text,
	"isolate" text,
	"background" text,
	"ra" numeric,
	"folder" text,
	"created" text,
	"modified" text,
	"note" text
);
-- Table 
CREATE TABLE IF NOT EXISTS "log" (
	"logId" integer PRIMARY KEY AUTOINCREMENT,
	"pipelineId" integer NOT NULL,
	"timeStamp" text,
	"logCode" integer,
	"logMsg" text,
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId")
);
-- Table 
CREATE TABLE IF NOT EXISTS "scaffold" (
	"pipelineId" integer NOT NULL,
	"scaffoldId" integer NOT NULL,
	"geneId" integer NOT NULL,
	"length" integer,
	"coverage" numeric,
	"geneStart" integer,
	"geneEnd" integer,
	PRIMARY KEY("pipelineId", "scaffoldId", "geneId")
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId")
);
-- Table 
CREATE TABLE IF NOT EXISTS "argInfo" (
	"pipelineId" integer NOT NULL,
	"qseqid" text NOT NULL,
	"geneId" integer NOT NULL,
	"sseqid" text,
	"qlen" integer,
	"slen" integer,
	"pident" numeric,
	"nident" integer,
	"length" integer,
	"evalue" numeric,
	"bitscore" numeric,
	"qstart" integer,
	"qend" integer,
	"gene" text,
	"subtype" text,
	"coverage" numeric,
	"type" text,
	"grp" integer,
	PRIMARY KEY("pipelineId", "geneId", "qseqId")
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId")
);
-- Table
CREATE TABLE IF NOT EXISTS "blastOut" (
	"pipelineId" integer NOT NULL,
	"qseqid" text NOT NULL,
	"sallacc" text,
	"staxids" integer,
	"sscinames" text,
	"salltitles" text,
	"qlen" integer,
	"slen" integer,
	"qstart" integer,
	"qend" integer,
	"sstart" integer,
	"send" integer,
	"bitscore" numeric,
	"length" integer,
	"pident" numeric,
	"nident" integer,
	"qcovs" integer,
	"rank" integer,
	"singleTop" integer,
	PRIMARY KEY("pipelineId", "qseqId", "sallacc")
	FOREIGN KEY("pipelineId") REFERENCES "pipeline"("pipelineId")
);
COMMIT;
