import json, os, stat, sys
from dataclasses import asdict
from pathlib import Path
import pytest
from tr_calling_pipeline.errors import CommandExecutionError, CommandTimeoutError, ExpectedOutputMissingError
from tr_calling_pipeline.execution import CommandSpec, OutputDeclaration, execute
from tr_calling_pipeline.provenance import canonical_digest, file_identity
from tr_calling_pipeline.stages import ResumeReason, resume_eligibility, select_stages
from tr_calling_pipeline.tools import Tool, ToolId, ToolStatus, detect_version, resolve_tool

def fake(tmp_path, body):
 p=tmp_path/'tool with spaces.py'; p.write_text('#!'+sys.executable+'\n'+body); p.chmod(p.stat().st_mode|stat.S_IXUSR); return p

def test_tool_resolution_and_version(tmp_path, monkeypatch):
 p=fake(tmp_path,"import sys; print('fake version 1.2.3', file=sys.stderr)\n")
 t=resolve_tool(Tool(ToolId.VAMOS,'VAMOS',str(p),True),tmp_path)
 assert t.configured_executable==str(p) and t.status is ToolStatus.AVAILABLE
 assert detect_version(t).detected_version=='1.2.3'
 missing=resolve_tool(Tool(ToolId.VAMOS,'VAMOS','absent',False),tmp_path,path='')
 assert missing.status is ToolStatus.MISSING_OPTIONAL

def test_execution_success_redaction_and_arguments(tmp_path):
 p=fake(tmp_path,"import json,os,sys; print(json.dumps(sys.argv[1:])); print(os.environ['TOKEN'],file=sys.stderr); open(sys.argv[1],'w').write('ok')\n")
 out=tmp_path/'out file'; spec=CommandSpec('c','s','VAMOS',(str(p),str(out),'a; $(bad)',''),str(tmp_path.resolve()),{'TOKEN':'secret'},('TOKEN',),declared_outputs=(OutputDeclaration('o',str(out)),))
 record=execute(spec,Tool(ToolId.VAMOS,'VAMOS',str(p),resolved_executable=str(p)),tmp_path/'record.json',tmp_path/'logs')
 assert record['status']=='SUCCEEDED' and record['argv'][-1]=='' and record['environment_overrides']['TOKEN']=='<REDACTED>'
 assert record['output_file_identities'][0]['sha256']==file_identity(out).sha256

def test_fail_timeout_and_missing(tmp_path):
 bad=fake(tmp_path,"import sys; sys.exit(7)\n"); tool=Tool(ToolId.VAMOS,'VAMOS',str(bad),resolved_executable=str(bad))
 with pytest.raises(CommandExecutionError): execute(CommandSpec('bad','s','VAMOS',(str(bad),),str(tmp_path.resolve())),tool,tmp_path/'bad.json',tmp_path/'logs')
 sleeper=fake(tmp_path,"import time; time.sleep(5)\n")
 with pytest.raises(CommandTimeoutError): execute(CommandSpec('slow','s','VAMOS',(str(sleeper),),str(tmp_path.resolve()),timeout_seconds=.1),tool,tmp_path/'slow.json',tmp_path/'logs')
 ok=fake(tmp_path,"pass\n")
 with pytest.raises(ExpectedOutputMissingError): execute(CommandSpec('missing','s','VAMOS',(str(ok),),str(tmp_path.resolve()),declared_outputs=(OutputDeclaration('x',str(tmp_path/'no')),)),tool,tmp_path/'missing.json',tmp_path/'logs')

def test_dry_run_and_stage_resume(tmp_path):
 p=tmp_path/'output'; p.write_text('x'); ident=asdict(file_identity(p)); tool=[{'tool_id':'VAMOS','resolved_executable':'/x','detected_version':'1','execution_mode':'NATIVE','container_digest':None}]
 prior={'record_schema_version':'1.0','status':'SUCCEEDED','configuration_digest':canonical_digest({'a':1}),'input_file_identities':[],'output_file_identities':[ident],'tool_identities':tool}
 assert resume_eligibility(prior,canonical_digest({'a':1}),[],[ident],tool) is ResumeReason.RESUME_ALLOWED
 p.write_text('y'); assert resume_eligibility(prior,canonical_digest({'a':1}),[],[asdict(file_identity(p))],tool) is ResumeReason.OUTPUT_CHANGED
 assert len(select_stages('02_align_assembly','05_run_straglr'))==4
 with pytest.raises(ValueError): select_stages('05_run_straglr','02_align_assembly')
