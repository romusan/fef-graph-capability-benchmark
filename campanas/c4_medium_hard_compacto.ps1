# Campana C4: medium-hard bajo geometria COMPACTA (180 mm).
# Replica exactamente la campana medium-hard publicada (paper_stats_20260521_131507)
# cambiando UNA sola variable: el perfil de empaquetamiento. Es el experimento
# que el propio paper declara que falta para aislar el efecto del envolvente.
param([int[]]$Semillas = (1..30), [int]$Max = 3)
$py   = "C:\Users\claudia.parra\.cadquery\cq311\python.exe"
$raiz = "C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\codigos"
$log  = Join-Path $raiz "campana_c4.log"
$est  = Join-Path $raiz "campana_c4_estado.json"
$salida = Join-Path $raiz "c4_mh_compacto"
New-Item -ItemType Directory -Force $salida | Out-Null
"[$(Get-Date -f 'HH:mm:ss')] C4: $($Semillas.Count) semillas, $Max en paralelo" |
  Out-File $log -Encoding utf8 -Append

$cola = New-Object System.Collections.Queue
foreach ($s in $Semillas) { $cola.Enqueue($s) }
$activos = @(); $hechas = 0; $fallos = 0; $inicio = Get-Date

while ($cola.Count -gt 0 -or $activos.Count -gt 0) {
  while ($activos.Count -lt $Max -and $cola.Count -gt 0) {
    $s = $cola.Dequeue()
    $j = Start-Job -ScriptBlock {
      param($py, $raiz, $salida, $seed)
      Set-Location $raiz
      # --- identico al medium-hard publicado ---
      $env:FEF_TARGET_MODE = "medium_hard"
      $env:FEF_TOPOLOGIES  = "wattI"
      $env:FEF_N_RANDOM = "900"; $env:FEF_SHAPE_SEEDS = "36"
      $env:FEF_TOPK_REFINE = "8"; $env:FEF_LOCAL_MAXITER = "100"
      $env:FEF_LS_MAX_NFEV = "800"
      $env:FEF_ELITE_NEIGHBORHOODS = "8"; $env:FEF_ELITE_NEIGHBOR_SAMPLES = "6"
      $env:FEF_RUN_DE = "1"; $env:FEF_DE_RESTARTS = "2"
      $env:FEF_DE_MAXITER = "55"; $env:FEF_DE_POPSIZE = "8"
      $env:FEF_RUN_LS_POLISH = "1"; $env:FEF_PHASE_STRIDE = "1"
      $env:FEF_RANDOM_SAMPLER = "lhs"
      $env:FEF_ACCEPT_RMS_MM = "1.0"; $env:FEF_ACCEPT_MAX_MM = "2.0"
      # --- la UNICA diferencia: el envolvente ---
      $env:FEF_GEOMETRY_PROFILE = "compact"
      $env:FEF_MAX_PRINTABLE_LINK_MM = "180.0"
      $env:FEF_SEED = "$seed"
      $env:FEF_OUTPUT_ROOT = $salida
      $env:FEF_OUT_PREFIX = "c4"
      $env:FEF_OUT_TAG = ("mh_compact_seed{0:d3}" -f $seed)
      $env:FEF_SAVE_FIGURES = "0"; $env:FEF_EXPORT_STL = "0"
      $env:FEF_EXPORT_STEP = "0"; $env:FEF_EXPORT_SVG = "0"
      $env:FEF_EXPORT_ZIP = "0"
      & $py .\run_single_synthesis.py 2>&1 | Out-String
    } -ArgumentList $py, $raiz, $salida, $s
    $j | Add-Member -NotePropertyName Etiqueta -NotePropertyValue "seed$s"
    $activos += $j
    "[$(Get-Date -f 'HH:mm:ss')] lanzada  $($j.Etiqueta)" | Add-Content $log -Encoding utf8
  }
  Start-Sleep -Seconds 15
  $vivos = @()
  foreach ($j in $activos) {
    if ($j.State -eq "Running") { $vivos += $j; continue }
    $sal = Receive-Job $j 2>&1 | Out-String
    if ($j.State -eq "Failed" -or $sal -match "Traceback") { $fallos++; $mk = "FALLO " }
    else { $hechas++; $mk = "ok    " }
    "[$(Get-Date -f 'HH:mm:ss')] $mk $($j.Etiqueta)" | Add-Content $log -Encoding utf8
    Remove-Job $j -Force
  }
  $activos = $vivos
  @{ hechas = $hechas; fallos = $fallos; activas = $activos.Count
     pendientes = $cola.Count
     minutos = [math]::Round(((Get-Date) - $inicio).TotalMinutes, 1) } |
    ConvertTo-Json | Out-File $est -Encoding utf8
}
"[$(Get-Date -f 'HH:mm:ss')] FIN: $hechas ok, $fallos fallos, $([math]::Round(((Get-Date)-$inicio).TotalMinutes,1)) min" |
  Add-Content $log -Encoding utf8
