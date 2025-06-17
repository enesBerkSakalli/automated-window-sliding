# Logs Directory

This directory contains all Nextflow execution logs and work files to keep the project root clean.

## Directory Structure

- `work/` - Nextflow work directory containing task execution files
- `trace/` - Execution reports and traces
  - `execution_timeline.html` - Timeline of pipeline execution
  - `execution_report.html` - Detailed execution report
  - `execution_trace.txt` - Raw trace data
  - `pipeline_dag.html` - Pipeline DAG visualization
- `process/` - Individual process logs (if needed)

## Cleanup

To clean up logs from previous runs:

```bash
rm -rf logs/work/*
rm -f logs/trace/*
```

## Note

These directories are ignored by git (.gitignore) to prevent committing large log files and temporary work data.
