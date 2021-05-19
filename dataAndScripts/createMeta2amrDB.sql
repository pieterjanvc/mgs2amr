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
-- Table that stores all ARG
CREATE TABLE IF NOT EXISTS "ARG" (
    "geneId" integer primary key,
	"refId" text,
	"gene" test,
	"subtype" text,
	"info" text,
	"nBases" integer,
	"clusterNr" integer
);
-- Table that stores all ARG detected in a sample
CREATE TABLE IF NOT EXISTS "detectedARG" (
	"pipelineId" integer NOT NULL,
	"runId" integer NOT NULL,
    "geneId" integer NOT NULL,
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
COMMIT;
