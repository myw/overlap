Results of running the calculation:

```python
{
  # Complement strand to the lambda sequences
  # Variant with 434 (coding strand) starting upstream
  '434': [
    # end of 434/start of lambda     end of lambda/start of next 434
    [(),                             ()],
    [(Seq('ATTGT', DNAAlphabet()),), (Seq('AC', DNAAlphabet()),)],
    [(),                             ()]
  ],

  # Variant with lambda complement strand starting upstream
  'lambda': [
    # end of lambda/start of 434   end of 434/start of next lambda
    [(),                           ()],
    [(Seq('AC', DNAAlphabet()),),  ()],
    [(),                           ()]
  ]
}
```

For the setup to be feasible, the first two entires in each variant should have enough
overlap so that the coding parts of the lambda strand should be no more than 8 bp
apart. Currently, there is not enough overlap to allow for this.
