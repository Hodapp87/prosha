use std::rc::Rc;

/*
struct R<'a> {
    b: &'a dyn Fn() -> R<'a>,
}

#[derive(Copy, Clone)]
struct Foo {}
impl<'a> Foo {
    // These are valid, but not especially useful (if I am
    // transferring ownership then I cannot have any branching):
    fn fn1(self) -> R<'a> {
        R { b: & move || self.fn1() }
    }
    fn fn2(self) -> R<'a> {
        R { b: &|| self.fn2() }
    }
}
*/

// Below (using box instead of a trait object) follows similar rules:
struct S<'a> {
    b: Box<dyn Fn() -> S<'a>>,
}
#[derive(Copy, Clone)]
struct Foo2 {}
impl<'a> Foo2 {
    fn fn1(self) -> S<'a> {
        S { b: Box::new(move || self.fn1()) }
    }
    // Not valid (error[E0373]: closure may outlive the current
    // function, but it borrows `self`, which is owned by the current
    // function):
    //fn fn2(self) -> S<'a> {
    //    S { b: Box::new(|| self.fn2()) }
    //}
    // Not valid:
    //fn fn3(&self) -> S<'a> {
    //    S { b: Box::new(move || self.fn3()) }
    //}
    // Not valid:
    //fn fn4(&self) -> S<'a> {
    //    S { b: Box::new(|| self.fn4()) }
    //}
}

struct T<'a> {
    b: Rc<dyn Fn() -> T<'a> + 'a>,
}
#[derive(Copy, Clone)]
struct Foo3 {}
impl<'a> Foo3 {
    fn fn1(self) -> T<'a> {
        T { b: Rc::new(move || self.fn1()) }
    }
    // Not valid (E0373):
    //fn fn2(self) -> T<'a> {
    //    T { b: Rc::new(|| self.fn2()) }
    //}
    // Not valid:
    //fn fn3(&self) -> T<'a> {
    //    T { b: Rc::new(move || self.fn3()) }
    //}
    // Not valid:
    //fn fn4(&self) -> T<'a> {
    //    T { b: Rc::new(|| self.fn4()) }
    //}
    // But this is now valid because T can be cloned:
    fn fn5(self) -> (T<'a>, T<'a>) {
        let p = Rc::new(move || self.fn1());
        let p2 = p.clone();
        (T { b: p }, T { b: p2 })
    }
}

// Further, this is now valid too (lifetimes removed):
struct U {
    b: Rc<dyn Fn() -> U>,
}
#[derive(Copy, Clone)]
struct Foo4 {}
impl Foo4 {
    fn fn1(self) -> U {
        U { b: Rc::new(move || self.fn1()) }
    }
    fn fn5(self) -> (U, U) {
        let p = Rc::new(move || self.fn1());
        let p2 = p.clone();
        (U { b: p }, U { b: p2 })
    }
}

// I can get rid of Copy/Clone if I use FnOnce:
struct V {
    b: Rc<dyn FnOnce() -> V>,
}
struct Foo5 {}
impl Foo5 {
    fn fn1(self) -> V {
        V { b: Rc::new(move || self.fn1()) }
    }
    fn fn2(self) -> (V, V) {
        let p = Rc::new(move || self.fn1());
        let p2 = p.clone();
        (V { b: p }, V { b: p2 })
    }
    // and then either kind is fine:
    fn fn3(self) -> V {
        V { b: Rc::new(|| self.fn3()) }
    }
    fn fn4(self) -> (V, V) {
        let p = Rc::new(|| self.fn3());
        let p2 = p.clone();
        (V { b: p }, V { b: p2 })
        // but this confuses me a bit. doesn't this then let me call
        // an FnOnce... more than once?
    }
}

