# Campana C2: amplia Bplus (QD+surrogate+refinement) y Estrict (linea base)
# de 3 a 15 semillas sobre el objetivo medio, geometria compacta.
# Las semillas 101/202/303 ya existen y NO se tocan.
$ErrorActionPreference = "Continue"
$py    = "C:\Users\claudia.parra\.cadquery\cq311\python.exe"
$raiz  = "C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\qd_experiments"
$log   = Join-Path $raiz "campana_c3.log"
$est   = Join-Path $raiz "campana_c3_estado.json"
$MAX   = 3
$semillas = 101,202,303,404,505,606,707,808,909,1010,1111,1212,1313,1414,1515

$tareas = @()
foreach ($s in $semillas) {
  $tareas += ,@("run_bplus_matched.py", $s)
}
$total = $tareas.Count
"[$(Get-Date -f 'HH:mm:ss')] campana C3: $total corridas, $MAX en paralelo" | Out-File $log -Encoding utf8

$cola = New-Object System.Collections.Queue
foreach ($t in $tareas) { $cola.Enqueue($t) }
$activos = @(); $hechas = 0; $fallos = 0; $inicio = Get-Date

while ($cola.Count -gt 0 -or $activos.Count -gt 0) {
  while ($activos.Count -lt $MAX -and $cola.Count -gt 0) {
    $t = $cola.Dequeue()
    $j = Start-Job -ScriptBlock {
      param($py, $raiz, $script, $seed)
      Set-Location $raiz
      & $py $script $seed 2>&1 | Out-String
    } -ArgumentList $py, $raiz, $t[0], $t[1]
    $j | Add-Member -NotePropertyName Etiqueta -NotePropertyValue "$($t[0]) seed$($t[1])"
    $activos += $j
    "[$(Get-Date -f 'HH:mm:ss')] lanzada  $($j.Etiqueta)" | Add-Content $log -Encoding utf8
  }
  Start-Sleep -Seconds 15
  $vivos = @()
  foreach ($j in $activos) {
    if ($j.State -eq "Running") { $vivos += $j; continue }
    $sal = Receive-Job $j 2>&1 | Out-String
    if ($j.State -eq "Failed" -or $sal -match "Traceback") { $fallos++; $mk = "FALLO  " }
    else { $hechas++; $mk = "ok     " }
    "[$(Get-Date -f 'HH:mm:ss')] $mk $($j.Etiqueta)  $((($sal -split "`n") | Select-Object -Last 2) -join ' ')" |
      Add-Content $log -Encoding utf8
    Remove-Job $j -Force
  }
  $activos = $vivos
  $min = [math]::Round(((Get-Date) - $inicio).TotalMinutes, 1)
  @{ total = $total; hechas = $hechas; fallos = $fallos
     activas = $activos.Count; pendientes = $cola.Count
     minutos = $min; pid_lanzador = $PID } |
    ConvertTo-Json | Out-File $est -Encoding utf8
}
"[$(Get-Date -f 'HH:mm:ss')] FIN: $hechas ok, $fallos fallos, $([math]::Round(((Get-Date)-$inicio).TotalMinutes,1)) min" |
  Add-Content $log -Encoding utf8
