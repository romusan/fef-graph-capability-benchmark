"""Second audit: exact revised-PDF line locations and evidence-specific replies."""
import json,re,hashlib
from pathlib import Path
HERE=Path(__file__).parent
KEYS={
'Introduction':'intro','Model identity and geometry':'models','Targets and geometric error':'targets',
'Cycle and transmission checks':'validation','Search and exact computational budget':'search',
'Controlled quality-diversity experiment':'qd','Statistics and reproducibility':'statistics',
'What the archived successes establish':'historical','Factorial target, profile and budget comparison':'factorial',
'Independent policy and additional paths':'independent','Recovery and discretization checks':'recovery',
'Penalty and transmission quality':'penalty','Quality-diversity controls':'qd_results',
'Discussion':'discussion','Conclusions':'conclusions','Complete effective parameter bounds':'bounds',
'Descriptive sampled point envelopes':'envelopes'}
MAP={
'1.1':['intro','independent','discussion'],'1.2':['search','factorial','bounds'],
'1.3':['models','independent','recovery'],'1.4':['models','envelopes'],'1.5':['targets','validation','recovery'],
'1.6':['statistics','factorial'],'1.7':['historical'],'1.8':['statistics'],'1.9':['qd','qd_results'],
'2.1':['intro','models','search'],'2.2':['models','search'],'2.3':['statistics','historical','factorial'],
'2.4':['models','targets','search'],'2.5':['intro','discussion'],'2.6':['ack'],
'3.1':['models','bounds'],'3.2':['validation','penalty'],'3.3':['targets'],'3.4':['targets','recovery'],
'3.5':['search','factorial'],'3.6':['historical'],'3.7':['statistics','factorial'],
'3.8':['qd','qd_results'],'3.9':['qd','qd_results'],'3.10':['statistics','qd_results']}

def instrument_body(body):
    matches=list(re.finditer(r'\\(?:sub)?section\{([^}]+)\}(?:\\label\{[^}]+\})?',body))
    out=body[:matches[0].start()]
    for i,m in enumerate(matches):
        end=matches[i+1].start() if i+1<len(matches) else len(body)
        segment=body[m.end():end]
        if m[1] not in KEYS:
            out+=m[0]+segment
            continue
        key=KEYS[m[1]]
        paragraphs=list(re.finditer(r'(?m)^[A-Za-z][^\n]*',segment))
        if not paragraphs:raise ValueError('No prose anchors for '+key)
        first=paragraphs[0].start();last=paragraphs[-1].end()
        segment=segment[:last]+r'\linelabel{rev:'+key+':end}'+segment[last:]
        segment=segment[:first]+r'\linelabel{rev:'+key+':start}'+segment[first:]
        out+=m[0]+segment
    return out

def strengthen_responses(responses):
    out=[]
    for ident,answer,changes,evidence in responses:
        changes=changes.replace('Results historical context','Results What the archived successes establish')
        if ident=='1.3':
            answer+=' We directly re-simulated the retained RC-G and B-O2 minima. Their reported errors are reproduced from real curves, not assigned on closure failure: invalid evaluations carry no RMS. The NACA0012 medians round to 3.929 mm because both searches favor almost rectilinear paths (median minor/major singular-value ratios 0.000992 and 0.000645). For CMA-ES/medium, 28/30 RC-G and 27/30 B-O2 paths have ratio below 0.1; the corresponding dense closure counts are 30/30 and 19/30. We provide median-nearest path plots, exact vectors and all 240 diagnostic records. B-O2/easy errors span 4.306--18.997 mm; the 6.213 median is not a fixed fill value. Own-curve recovery tests a different feasible target and does not establish coverage of the other targets.'
            changes+=' Results Independent policy and additional paths; diagnostic figure.'
            evidence+=' 15_auditoria_reenvio/model_diagnostics.json; selected_paths_diagnosis.png; example coordinate CSVs.'
        if ident=='1.4':
            answer+=' We additionally report an actual geometric descriptor: the dense sampled axis-aligned envelope of all represented moving points, fixed anchors and tracer, and the maximum radius from the crank ground anchor. The new appendix covers all 960 factorial minimum-RMS selections at the base budget, with dense-closure failures explicitly counted. This is a post-hoc point-envelope measurement, not an imposed bound or a collision certificate; body thickness is not modeled. We therefore retain the term parameter-bound profile rather than interpreting this measurement as packaging feasibility.'
            changes+=' Appendix Descriptive sampled point envelopes.'
            evidence+=' 15_auditoria_reenvio/envelope_diagnostics.json and envelope_diagnostics_summary.csv.'
        if ident in ['1.7','3.6']:
            answer+=' An exploratory Pearson 3-by-2 homogeneity diagnostic gives chi-square 2.219, df=2, p=0.3297 (minimum expected count 7.67). This does not demonstrate identical batch rates or zero variability.'
            evidence+=' 15_auditoria_reenvio/model_diagnostics.json, batch_heterogeneity.'
        if ident=='2.5':
            answer+=' We now use the peer-reviewed Plecnik and McCarthy (2016) article directly for the Watt/Stephenson distinction and remove the teaching-slide citation. No graph-search performance comparison was run: graph synthesis changes the structural search space, whereas the present experiment holds eight specified models fixed. The contribution is a controlled validation and recoverability study, with no claim of outperforming structural synthesis.'
        if ident=='1.2':answer+=' A campaign inventory table now explicitly sums 960 factorial, 480 independent-policy, 1200 additional/easy/reference-target and 240 penalty-ablation searches to 2880, separately from repetitions, recovery and QD.'
        out.append((ident,answer,changes,evidence))
    return out

def location_text(ident,paper):
    path=HERE/'results/revised_line_locations.json'
    if not path.exists():return 'Locations will be regenerated from the compiled revised manuscript.'
    d=json.loads(path.read_text())
    if d['source_sha256']!=hashlib.sha256(paper.read_bytes()).hexdigest():return 'Locations will be regenerated from the compiled revised manuscript.'
    parts=[]
    for key in MAP[ident]:
        r=d['locations'][key]
        parts.append(f"pp. {r['first_page']}--{r['last_page']}, lines {r['first_line']}--{r['last_line']}")
    return '; '.join(parts)+'.'
