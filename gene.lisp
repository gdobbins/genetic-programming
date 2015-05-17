(in-package :cl-user)
(eval-when (:compile-toplevel :load-toplevel :execute)
  ;(ql:quickload :series)
  (ql:quickload :alexandria)
  (ql:quickload :lparallel))
(defpackage :genetic-programming
  (:use :common-lisp :alexandria :lparallel)
  (:export :gene-prog :make-population :population-genes :population-gene-fitness-list :population-gene-adjusted-fitness-list :population-number-of-genes :population-average-fitness :population-maximum-fitness :population-minimum-fitness :population-best-gene-so-far :population-best-gene-so-far-fitness :make-functions :functions-function-list :functions-arg-num-list :simulated-annealing))
(in-package :genetic-programming)

;(defclass gene ()
;  ((gene :accessor gene :initarg :gene)
;   (raw-fitness :accessor raw-fitness)
;   (adjusted-fitness :accessor adjusted-fitness)))

(defvar *number-of-threads* 4)
(setf lparallel:*kernel* (lparallel:make-kernel *number-of-threads*))



(defstruct population
  (genes nil :type (or cons null))
  (gene-fitness-list nil :type (or cons null))
  (gene-adjusted-fitness-list nil :type (or cons null))
  (number-of-genes 0 :type fixnum)
  (average-fitness 0 :type number)
  (maximum-fitness 0 :type number)
  (minimum-fitness 0 :type number)
  (best-gene-so-far nil :type (or cons null))
  (best-gene-so-far-fitness 0 :type number))

(defstruct functions
  (function-list nil :type (or cons null))
  (arg-num-list nil :type (or cons null))
  (base-gene-copy-number 1 :type fixnum))

;(defmacro concurrent-lists-do (lists &body body)
;  (let* ((syms (loop for i in lists collect (gensym))) (len (length syms)))
;    `(do ,(loop for i from 0 to (1- len) collect `(,(nth i syms) ,(nth i lists) (rest ,(nth i syms))))
;	 ((or ,@(loop for i from 0 to (1- len) collect `(null ,(nth i syms)))))
;       (let ,(loop for i from 0 to (1- len) collect `(,(nth i lists) ,(nth i syms)))
;       ,@body))))

(defun tree-replace (tree find max-depth min-depth replace &rest vars)
  (declare (type cons tree) (type fixnum max-depth min-depth))
  (loop for i in tree collect
       (if (and (consp i) (<= 0 max-depth) (not #1=(find-if (lambda (x) (equalp i x)) (last vars))))
	   (apply 'tree-replace i find (1- max-depth) (1- min-depth) replace vars)
	   (if (por (equal i find) (and (<= 0 min-depth) #1#) (and (>= 0 max-depth) (consp i) (not #1#)))
	       (apply replace vars)
	       i))))

(defun tree-find-not (tree find)
  (declare (type cons tree))
  (let ((rep t))
    (declare (type boolean rep))
    (loop for i in tree do
	 (if rep
	     (if (consp i)
		 (setf rep (tree-find-not i find))
		 (if (equal i find)
		     (setf rep nil)))))
    rep))

(defun tree-depth (tree &optional (depth 0) (max-depth 0) (ret nil))
  (declare (type (or cons null) tree ret) (type fixnum depth max-depth))
  (if tree
      (if (consp (car tree))
	  (tree-depth (cdar tree) (1+ depth) (max max-depth (1+ depth)) (cons (list depth (cdr tree)) ret))
	  (tree-depth (cdr tree) depth max-depth ret))
      (if ret
	  (tree-depth (second (first ret)) (first (first ret)) max-depth (rest ret))
	  max-depth)))

(defun tree-atom-size (tree &optional (functions nil) (atoms 0) (ret nil))
  "Counts the number of atoms in a tree. Does not count () expressions as atoms, skips the function slot of the list, i.e. the first position of every list if functions is nil. Subtract one if gene is not wrapped in a list."
  (declare (type (or cons null) tree ret) (type boolean functions) (type fixnum atoms))
  (if tree
      (if (consp (car tree))
	  (tree-atom-size (cdar tree) functions (if functions (1+ atoms) atoms) (cons (cdr tree) ret))
	  (tree-atom-size (cdr tree) functions (1+ atoms) ret))
      (if ret
	  (tree-atom-size (first ret) functions atoms (rest ret))
	  atoms)))

(defun node-list (gene &optional (acc nil) (ret nil))
  (declare (type (or cons null) gene acc ret))
  (if gene
      (if (consp (car gene))
	  (node-list (cdar gene) (cons gene acc) (cons (cdr gene) ret))
	  (node-list (cdr gene) (cons gene acc) ret))
      (if ret
	  (node-list (first ret) acc (rest ret))
	  acc)))

(defgeneric create-random-gene (functions terminals &optional max-depth min-depth))

(defmethod create-random-gene ((functions functions) (terminals cons) &optional (max-depth 6) (min-depth 2))
  (declare (type fixnum max-depth min-depth))
  (let* ((len-fun (length #1=(functions-function-list functions))) (len-term (length terminals)) (len (+ len-fun len-term)) gene)
    (declare (type fixnum len-fun len-term len) (type (or null cons)))
    (flet ((create-new-sexp (num functions terminals)
	     (declare (type fixnum num) (type functions functions) (type cons terminals))
		 (if (< num len-fun)
		     (append (list (nth num #1#)) (make-list (nth num (functions-arg-num-list functions)) :initial-element :sexp))
		     (nth (- num len-fun) terminals))))
      (setf gene (create-new-sexp (random len-fun) functions terminals))
      (do () ((tree-find-not gene :sexp) gene)
	(setf gene (tree-replace gene :sexp max-depth min-depth (lambda (x y z) (create-new-sexp (random x) y z)) len functions terminals))))))

(defun create-random-full-gene (functions terminals &optional (max-depth 6) (min-depth 2))
  (declare (type functions functions) (type cons terminals) (type fixnum max-depth min-depth))
  (loop for i from 1 to (functions-base-gene-copy-number functions) collect
       (create-random-gene functions terminals max-depth min-depth)))

(defun initialize-population (current-population functions terminals generation)
  (declare (type population current-population) (type cons terminals) (type functions functions) (type fixnum generation))
  (assert (= 0 generation))
  (if (< (length #1=(population-genes current-population)) #2=(population-number-of-genes current-population))
      (tagbody (setf (population-genes current-population) (append (list (create-random-full-gene functions terminals)) (population-genes current-population))) (initialize-population current-population functions terminals generation))
;      (initialize-population (make-population :genes (append (list (create-random-gene functions terminals)) #1#) :number-of-genes #2#) functions terminals generation)
      current-population))

(defun adjust-to-linear-fitness (fitness-list max-fit min-fit avg-fit fit-list-length &optional (max-fit-multiplier 2))
  (declare (type cons fitness-list) (type number max-fit min-fit avg-fit) (type fixnum fit-list-length max-fit-multiplier))
  (if (= max-fit avg-fit) (make-list fit-list-length :initial-element 1)
      (let* ((a (/ (- max-fit-multiplier 1) (- max-fit avg-fit))) (b (* -1 (1+ (* a avg-fit)))) (result #1=(pmapcar (lambda (x) (declare (type number x)) (+ b (* a x))) :size fit-list-length fitness-list)))
	(declare (type number a b) (type cons result))
;	(setf a (/ (- max-fit-multiplier 1) (- max-fit avg-fit)) b (* -1 (1+ (* a avg-fit))))
;	#1=(setf result (pmapcar (lambda (x) (declare (type number x)) (+ b (* a x))) fitness-list)
	(if (psome #'minusp :size fit-list-length result)
	    (if (= avg-fit min-fit)
		(setf result (make-list fit-list-length :initial-element 1))
		(setf a (/ 1 (- avg-fit min-fit)) b (* -1 a min-fit) result #1#)))
	(pmap-into result (lambda (x) (declare (type number x)) (coerce x 'single-float)) :size fit-list-length result))))

(defun update-population-fitness (current-population fitness-function)
  (declare (type function fitness-function) (type population current-population))
  (labels ((upf-best-gene-helper (gene-list fitness-list best-gene-fitness)
	     (declare (type cons gene-list fitness-list) (type number best-gene-fitness))
	     (if (= (car fitness-list) best-gene-fitness)
		 (values (car gene-list) (car fitness-list))
		 (upf-best-gene-helper (cdr gene-list) (cdr fitness-list) best-gene-fitness))))
    (pmap-into #1=(population-gene-fitness-list current-population) fitness-function :size (population-number-of-genes current-population) (population-genes current-population))
    (setf (population-average-fitness current-population) (/ (preduce #'+ #1#) (population-number-of-genes current-population)))
    (setf (population-minimum-fitness current-population) (preduce #'min #1#))
    (setf (population-maximum-fitness current-population) (preduce #'max #1#))
    (setf (population-gene-adjusted-fitness-list current-population) (adjust-to-linear-fitness (population-gene-fitness-list current-population) (population-maximum-fitness current-population) (population-minimum-fitness current-population) (population-average-fitness current-population) (population-number-of-genes current-population)))
    (if (> (population-maximum-fitness current-population) (population-best-gene-so-far-fitness current-population))
	(multiple-value-bind (gene fitness) (upf-best-gene-helper (population-genes current-population) (population-gene-fitness-list current-population) (population-maximum-fitness current-population))
	  (setf (population-best-gene-so-far current-population) gene (population-best-gene-so-far-fitness current-population) fitness)))))

(defun cmp-helper (gene-list fitness-list selector-number)
  (declare (type single-float selector-number) (type cons gene-list fitness-list) (optimize (speed 3) (safety 0) (compilation-speed 0) (debug 0)))
  (if (> (the single-float (car fitness-list)) selector-number)
      (car gene-list)
      (if (cdr gene-list)
	  (cmp-helper (cdr gene-list) (cdr fitness-list) (- selector-number (the single-float (car fitness-list))))
	  (cmp-helper gene-list fitness-list (1- selector-number)))))

(defun create-mating-pool (current-population)
  (declare (type population current-population))
  (plet ((total-fitness (reduce (lambda (x y) (declare (type single-float x y) (optimize (speed 3) (safety 0))) (handler-case (+ x y) (floating-point-overflow () most-positive-single-float))) (population-gene-adjusted-fitness-list current-population))))
    (declare (type single-float total-fitness))
    (make-population :number-of-genes (population-number-of-genes current-population) :best-gene-so-far (population-best-gene-so-far current-population) :best-gene-so-far-fitness (population-best-gene-so-far-fitness current-population)
		     :average-fitness (population-average-fitness current-population)
		     :maximum-fitness (population-maximum-fitness current-population)
		     :minimum-fitness (population-minimum-fitness current-population)
		     :gene-fitness-list (population-gene-fitness-list current-population)
		     :gene-adjusted-fitness-list (population-gene-adjusted-fitness-list current-population)
		     :genes (pmapcar (lambda (x) (declare (ignore x)) (copy-tree (cmp-helper (population-genes current-population) (population-gene-adjusted-fitness-list current-population) (random total-fitness)))) (population-genes current-population)))))

;		     (loop for i from 1 to (population-number-of-genes current-population) collect (copy-tree (cmp-helper (population-genes current-population) (population-gene-adjusted-fitness-list current-population) (random total-fitness)))))))
(defun bmp-helper (gene-list functions terminals pairs-remaining mutate-percentage &optional (max-depth 50))
  (declare (type fixnum pairs-remaining max-depth) (type cons gene-list terminals) (type single-float mutate-percentage)
	   (type functions functions))
  (if (<= 2 pairs-remaining)
      (let ((base-gene-number (random (functions-base-gene-copy-number functions))))
	(plet ((sub-one (random-elt (node-list (list (nth base-gene-number (first gene-list)))))) (sub-two (random-elt (node-list (list (nth base-gene-number (second gene-list)))))))
	  (declare (type cons sub-one sub-two))
	  (rotatef (car sub-one) (car sub-two))
	  (plet ((tree-depth-of-second-gene-list (tree-depth (second gene-list))))
	    (if (< max-depth (the fixnum (tree-depth (first gene-list))))
		(if #1=(< max-depth (the fixnum tree-depth-of-second-gene-list))
		    (rotatef (car sub-one) (car sub-two))
		    (setf (car sub-one) (car sub-two)))
		(if #1#
		    (setf (car sub-two) (car sub-one)))))
	  (bmp-helper (cddr gene-list) functions terminals (- pairs-remaining 2) mutate-percentage max-depth)))
      (let ((n (length gene-list)))
	(declare (type fixnum n))
	(dotimes (i (the fixnum (round (+ (* n mutate-percentage (- 1 mutate-percentage) (gaussian-random)) (* n mutate-percentage)))))
	  (plet ((sub-gene (random-elt gene-list)))
	    (plet ((sub-one (random-elt (node-list sub-gene))) (sub-two (list (create-random-gene functions terminals -1 1))))
	      (rotatef (car sub-one) (car sub-two))
	      (if (< max-depth (the fixnum (tree-depth sub-gene)))
		  (rotatef (car sub-one) (car sub-two)))))))))

(defun breed-mating-pool (mating-pool functions terminals &optional (breed-percentage 0.9) (mutate-percentage .01))
  (declare (type population mating-pool) (type functions functions) (type cons terminals) (type single-float breed-percentage mutate-percentage))
  (bmp-helper (population-genes mating-pool) functions terminals (round (* (population-number-of-genes mating-pool) breed-percentage)) mutate-percentage)
    mating-pool)

(declaim (inline update-history-report))
(defun update-history-report (current-population generation print-freq)
  (declare (type fixnum generation print-freq))
  (if (zerop (mod generation print-freq))
      (format t "~&--------------------------------------------------~%Generation: ~d     Stardate: ~{~d-~d-~d  ~d:~2,'0d:~2,'0d~}~%Maximum Fitness: ~d   Minimum Fitness: ~d   Average Fitness: ~d   Best Fitness: ~d   "
	      generation (nreverse (subseq (multiple-value-list (get-decoded-time)) 0 6)) (population-maximum-fitness current-population) (population-minimum-fitness current-population) (population-average-fitness current-population) (population-best-gene-so-far-fitness current-population))))

(let ((current-population (make-population)) (mating-pool (make-population))
      (terminals nil) (functions (make-functions)) (print-freq 1)
      (fitness-function #'(lambda (x) (eval (car x))))
      (generation 0))
  (declare (type fixnum generation print-freq) (type (or null cons) terminals) (type function fitness-function) (type functions functions) (type population current-population mating-pool))
  (defun gene-prog (&key (operation :next-generation) (generations 50);for progressing to the next evolution
		      (delete-duplicates nil) (population-size 100) (init-functions functions) (init-terminals terminals) (fit-function fitness-function) (print-frequency print-freq);for during initialization, whether to check for duplicates
		      (file nil)) ;for saving and loading
    (declare (type fixnum generations population-size print-frequency) (type keyword operation) (type boolean delete-duplicates) (type functions init-functions) (type cons init-terminals) (type function fit-function) (type (or null (simple-array character (*)))))
    (case operation
      (:next-generation
       (update-population-fitness current-population fit-function)
       (setf mating-pool (create-mating-pool current-population))
       (setf current-population (breed-mating-pool mating-pool functions terminals))
       (update-history-report current-population generation print-freq)
       (setf generation (1+ generation)))
      (:best-of-generation
       (labels ((bog-helper (gene-list fitness-list max-fitness)
		  (if (= max-fitness (car fitness-list))
		      (car gene-list)
		      (bog-helper (cdr gene-list) (cdr fitness-list) max-fitness))))
	 (update-population-fitness current-population fit-function)
	 (values (bog-helper (population-genes current-population) (population-gene-fitness-list current-population) (population-maximum-fitness current-population)) (population-maximum-fitness current-population))))
      (:best-of-generation-continue
       (setf mating-pool (create-mating-pool current-population))
       (setf current-population (breed-mating-pool mating-pool functions terminals))
       (update-history-report current-population generation print-freq)
       (setf generation (1+ generation)))
      (:next-generations
       (dotimes (i generations)
	 (update-population-fitness current-population fitness-function)
	 (setf mating-pool (create-mating-pool current-population))
	 (setf current-population (breed-mating-pool mating-pool functions terminals))
	 (update-history-report current-population generation print-freq)
	 (setf generation (1+ generation))))
       (:initialize
       (setf terminals init-terminals functions init-functions
	     (population-number-of-genes current-population) population-size
	     (population-number-of-genes mating-pool) population-size)
       (setf fitness-function fit-function)
       (setf print-freq print-frequency)
       (initialize-population current-population functions terminals generation)
       (setf (population-gene-fitness-list current-population) (make-list (population-number-of-genes current-population) :initial-element 0))
       (if delete-duplicates (funcall 'gene-prog :operation :delete-initial-duplicates
				      :population-size population-size
				      :init-functions init-functions
				      :init-terminals init-terminals
				      :delete-duplicates t)))
      (:enlarge-population
       (setf (population-number-of-genes current-population) population-size
	     (population-number-of-genes mating-pool) population-size)
       (initialize-population current-population functions terminals 0)
       (setf (population-gene-fitness-list current-population) (make-list (population-number-of-genes current-population) :initial-element 0)))
      (:delete-initial-duplicates
       (let ((old-population-size (length (population-genes current-population))))
       (setf #1=(population-genes current-population) (delete-duplicates #1# :test #'tree-equal))
       (funcall 'gene-prog :operation :initialize
		:population-size old-population-size
		:init-functions init-functions
		:init-terminals init-terminals
		:print-frequency print-frequency
		:delete-duplicates (if (= old-population-size (length #1#)) nil t))))
      (:save (with-open-file (out file :direction :output :if-exists :supersede)
	       (prin1 current-population out)
	       (prin1 mating-pool out)
	       (prin1 terminals out)
	       (prin1 functions out)
	       (prin1 (list generation) out)
	       (prin1 (list print-freq) out)))
      (:load (with-open-file (in file :direction :input)
	       (setf current-population (read in))
	       (setf mating-pool (read in))
	       (setf terminals (read in))
	       (setf functions (read in))
	       (setf generation (car (read in)))
	       (setf print-freq (car (read in)))))
      (:change-terminals (setf terminals init-terminals))
      (:set-fit-function (setf fitness-function fit-function))
      (:mating-pool mating-pool)
      (:fitness-function fitness-function)
      (:terminals terminals)
      (:functions functions)
      (:generation generation)
      (:print-freq print-freq)
      (:current-population current-population))))

(defun boltzmann (e-i e-i+1 temp &optional (k 1.0d0))
  (declare (type real e-i e-i+1) (type double-float temp k))
  (exp (the double-float (/ (abs (- e-i e-i+1)) -1 k temp))))

(defun simulated-annealing (initial-value initial-temperature energy-function step-function &optional (min-t 2.0d-6) (mu-T 1.003d0) (k 1.0d0))
  (declare (type function energy-function step-function) (type double-float initial-temperature min-t mu-t k))
  (let ((value initial-value) (value-energy (funcall energy-function initial-value)))
    (do* ((candidate-value (funcall step-function value) (funcall step-function value)) (candidate-value-energy (funcall energy-function candidate-value) (funcall energy-function candidate-value)) (temp initial-temperature (/ temp mu-T))) ((< temp min-t) (if (< candidate-value-energy value-energy) (values candidate-value candidate-value-energy) (values value value-energy)))
      (declare (type double-float temp))
      (if (or (< candidate-value-energy value-energy) (< (random 1.0d0) (the double-float (boltzmann value-energy candidate-value-energy temp k))))
	  (setf value candidate-value value-energy candidate-value-energy)))))
	      
