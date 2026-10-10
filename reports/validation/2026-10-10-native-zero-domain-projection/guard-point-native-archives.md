The finalized guard-point two-loop fixed-dimension, two-loop Laurent, and failed three-loop native corpora were released by the root and acceptance owners before archival. The new `guard-point-native-archives.json` binds each original path, SHA256 and byte count to its deterministic gzip payload. The massive regression corpus is excluded.

All native binary programs, native JSON payloads of at least 8192 bytes, and the three sibling `*-resources.log` files are compressed. Predictions, configurations, inputs, resource summaries, comparisons, bindings and small native diagnostics remain readable at their original paths. Every original payload was hashed, compressed with gzip level 9 and timestamp zero, decompressed and verified before its raw copy was removed. Existing frozen evidence manifests are unchanged.

Run `python verify-guard-point-native-archives.py` from any directory to verify all archived and preserved bytes. Historical probe commands still name original raw paths. Restore a payload explicitly before running such a command:

```sh
gzip -dc ARCHIVE_PATH > ORIGINAL_PATH
```

Use the two exact paths from the corresponding manifest record. Restoring the raw copy does not change the archive or original hash. Archival records storage changes only; it does not turn the failed three-loop preparation into an accepted amplitude.
